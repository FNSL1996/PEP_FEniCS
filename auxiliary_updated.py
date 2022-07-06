#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 14 15:01:20 2021

@author: felipe
"""

'''
Auxiliary Functions: Modified to use mechanical convention
'''

import numpy as np
import dolfin


# Function sigma
# vector, scalar, scalar, scalar, scalar -> tensor CT(2)
def sigma(us,p,alpha,K,G):
    sigma_val = sigma_e(us,K,G) - alpha*p*dolfin.Identity(np.size(us))
    return sigma_val

# Stress effective
# vector, scalar, scalar -> tensor CT(2)
def sigma_e(us,K,G):
    sigma_e_val = 2.0*G*epsilon(us) + (K-2.0/3.0*G)*epsilon_v(us)*\
    dolfin.Identity(np.size(us))
    return sigma_e_val

# Tensor de Strain (pequeñas deformaciones)
# vector -> tensor CT(2)
def epsilon(us):
    epsilon_val = dolfin.sym(dolfin.grad(us))
    return epsilon_val

# Volumetric Strain
# vector -> scalar
def epsilon_v(us):
    epsilon_v_val = dolfin.div(us)
    return epsilon_v_val

# Volumetric Stress
def sigma_v(us,p,alpha,K,G):
    sigma_vol = dolfin.tr(sigma(us,p,alpha,K,G))/3.0
    return sigma_vol

# Stress deviatorico
def sigma_dev(us,p,alpha,K,G):
    sigma_dev = sigma(us,p,alpha,K,G) - sigma_v(us,p,alpha,K,G)
    return sigma_dev
    
   
# vector -> tensor CT(2)
def as_3D_tensor(X): 
    return dolfin.as_tensor([[X[0], X[3], X[4]],
                         [X[3], X[1], X[5]],
                          [X[4], X[5], X[2]]])
                          
def tensor2vector(X):
    return dolfin.as_vector([X[0,0], X[1,1], X[2,2], X[0,1], X[0,2], X[1,2]])
                          
# Isotropic hardening for yield strength
def sigma_y(sigy_0,eplas,H):
    return sigy_0 + H*eplas
                 
# Correction mapping for sigma
def sig_correction(u_trial,sig_prev,eplas_prev,K,G,sigy_0,H):
    # Parameters at first iteration   

    sig_trial = sig_prev + sigma_e(u_trial,K,G)
    s_trial = dolfin.dev(sig_trial)
    sig_eq_trial = dolfin.sqrt((3.0/2.0)*dolfin.inner(s_trial,s_trial)) # Equivalent von Mises stress
    sigy_trial = sigma_y(sigy_0,eplas_prev,H)
    f_trial = sig_eq_trial- sigy_trial
    # Plastic condition
    conditioner = dolfin.conditional(dolfin.le(f_trial,0),0,1)
    delta_lam = conditioner*f_trial/(3*G + H) # Plastic multiplier
    beta = 3*G*delta_lam/sig_eq_trial #Reduccion de respuesta deviatorica (buscar nombre en libro)
    N = conditioner*(s_trial)/sig_eq_trial
    sig_vol = sig_trial - s_trial
    s_now = (1 - beta)*s_trial
    sig_now = sig_vol + s_now
    e_now = (1/(2*G))*s_now + (epsilon_v(u_trial)/3)*dolfin.Identity(3)
    evol_now = dolfin.tr(e_now)
    edev_now = dolfin.dev(e_now)
    sig_eq_now = dolfin.sqrt((3.0/2.0)*dolfin.inner(s_now,s_now)) # Equivalent von Mises stress
    p_now = dolfin.tr(sig_now) # Mean stress    
    return sig_now, N, \
           beta, delta_lam, sig_eq_now, conditioner, p_now, evol_now, edev_now
           
# Tangent stress
def sigma_tang(K,G,H,u_now,beta,N):
    e_now = epsilon(u_now)
    return sigma_e(u_now,K,G) - 3*G*((3*G/(3*G + H)) - beta)*dolfin.inner(N,e_now)*(N) - beta*2*G*dolfin.dev(e_now) 


# Homogenization for DirichletBCS
def list_homogenize(BCS):
    return BCS.homogenize()
    
# Projection on to quadrature spaces: vonMises_plasticity.py
def local_project(v, V, dxm, u=None):
    dv = dolfin.TrialFunction(V)
    v_ = dolfin.TestFunction(V)
    a_proj = dolfin.inner(dv, v_)*dxm
    b_proj = dolfin.inner(v, v_)*dxm
    solver = dolfin.LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u is None:
        u = dolfin.Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return

# Combine vector for multiple materials
def combine_results(x,y):
    combined = np.zeros(len(x))
    for i in range(len(x)):
        if np.isnan(x[i]) == True:
            combined[i] = y[i]
        else:
            combined[i] = x[i]
    return combined

# Combine fields for two materials (can be modified for more)
# Domains are the integration domain
def combine_fields_original(fields, names, vectorspace, mesh, volumes, domains):
    output_fields = []
    for i in range(0,len(fields[0,:])):
        dofmap = vectorspace[i].dofmap()
        field1 = fields[0,i]
        field2 = fields[1,i]
        test_field = field1
        for cell in dolfin.cells(mesh):
            if volumes[cell.index()] != domains[0].subdomain_id():
                dofs = dofmap.cell_dofs(cell.index())
                for dof in dofs:
                  test_field.vector()[dof] = field2.vector()[dof]        
        test_field.rename(names[i],names[i])
        output_fields.append(test_field)
    return output_fields
        
# Combine fields for two materials (can be modified for more)
# Domains are the integration domain
def combine_field(fields, names, vectorspace, mesh, volumes, domains):
    output_fields = []
    dofmap = vectorspace.dofmap()
    field1 = fields[0]
    field2 = fields[1]
    test_field = field1
    for cell in dolfin.cells(mesh):
        if volumes[cell.index()] != domains[0].subdomain_id():
            dofs = dofmap.cell_dofs(cell.index())
            for dof in dofs:
              test_field.vector()[dof] = field2.vector()[dof]        
    test_field.rename(names[0],names[0])
    output_fields.append(test_field)
    return output_fields
        

def EigenSolver_DG(mesh,sigma,tensorspace,newarrayeig,vecarray,vectorspace,funarray,funspace):
    eigvalues = np.zeros((mesh.num_cells(),9))
    sig_aux = sigma.vector()[:]
    sig_eig = np.zeros(mesh.num_cells()*9)
    vec1 = np.zeros(mesh.num_cells()*3)
    vec2 = np.zeros(mesh.num_cells()*3)
    vec3 = np.zeros(mesh.num_cells()*3)
    phi = np.zeros(mesh.num_cells())
    for cell in dolfin.cells(mesh):
        dofmap = tensorspace.dofmap()
        cell_dofmap = np.array(dofmap.cell_dofs(cell.index()))
        sigma_cell = sig_aux[list(map(int,cell_dofmap))]
        [eig_cell,eigvec_cell] = np.linalg.eig(sigma_cell.reshape(3,3))
        idx = np.argsort(eig_cell) # Ascending
        eig_cell = eig_cell[idx]
        eigvec_cell = eigvec_cell[:,idx]
        sig_eig[list(map(int,cell_dofmap))[0::4]] = eig_cell
        # Eigenvector listing
        vec_dofmap = vectorspace.dofmap()
        veccell_dofmap = np.array(vec_dofmap.cell_dofs(cell.index()))
        vec1[list(map(int,veccell_dofmap))] = eigvec_cell[:,0]
        vec2[list(map(int,veccell_dofmap))] = eigvec_cell[:,1]
        vec3[list(map(int,veccell_dofmap))] = eigvec_cell[:,2]
        # Tensor shape
        fun_dofmap = funspace.dofmap()
        funcell_dofmap = np.array(fun_dofmap.cell_dofs(cell.index()))
        phi[list(map(int,funcell_dofmap))] = (eig_cell[1] - eig_cell[2])/(eig_cell[0] - eig_cell[2])       
    newarrayeig.vector().set_local(sig_eig)
    vecarray[0].vector().set_local(vec1)
    vecarray[1].vector().set_local(vec2)
    vecarray[2].vector().set_local(vec3)
    funarray.vector().set_local(phi)
    return 

def EigenSolver_DG_Order1(mesh,sigma,newarrayeig,vecarray,funarray):
    eigvalues = np.zeros(mesh.num_cells()*4*3)
    sig_aux = sigma.vector()[:]
    vec1 = np.zeros(mesh.num_cells()*4*3)
    vec2 = np.zeros(mesh.num_cells()*4*3)
    vec3 = np.zeros(mesh.num_cells()*4*3)
    phi = np.zeros(mesh.num_cells()*4)
    # Eigenvalues and eigenvectors
    [eig_ver,eigvec_ver] = np.linalg.eig(sig_aux.reshape([mesh.num_cells()*4,3,3]))
    eig_index = np.argsort(eig_ver, axis = 1)
    # Sorting eigenvalues
    eig_ver = np.take_along_axis(eig_ver,eig_index,axis=1)
    for i in range(0,mesh.num_cells()*4):
        phi[i] = (eig_ver[i,1] - eig_ver[i,2])/(eig_ver[i,0] - eig_ver[i,2])
        eigvec_range = eigvec_ver[i,:,:]
        eigvec_range = eigvec_range[:,eig_index[i,:]]
        vec1[3*i:3*i+3] = eigvec_range[:,0]
        vec2[3*i:3*i+3] = eigvec_range[:,1]
        vec3[3*i:3*i+3] = eigvec_range[:,2]
        eigvalues[3*i:3*i+3] = eig_ver[i,:]
    newarrayeig.vector().set_local(eigvalues)
    vecarray[0].vector().set_local(vec1)
    vecarray[1].vector().set_local(vec2)
    vecarray[2].vector().set_local(vec3)
    funarray.vector().set_local(phi)
    return 

def EigenSolver_CG(mesh,sigma,newarrayeig,vecarray,funarray):
    eigvalues = np.zeros(mesh.num_vertices()*3)
    sig_aux = sigma.vector()[:]
    vec1 = np.zeros(mesh.num_vertices()*3)
    vec2 = np.zeros(mesh.num_vertices()*3)
    vec3 = np.zeros(mesh.num_vertices()*3)
    phi = np.zeros(mesh.num_vertices())
    # Eigenvalues and eigenvectors
    [eig_ver,eigvec_ver] = np.linalg.eig(sig_aux.reshape([mesh.num_vertices(),3,3]))
    eig_index = np.argsort(eig_ver, axis = 1)
    # Sorting eigenvalues
    eig_ver = np.take_along_axis(eig_ver,eig_index,axis=1)
    for i in range(0,mesh.num_vertices()):
        phi[i] = (eig_ver[i,1] - eig_ver[i,2])/(eig_ver[i,0] - eig_ver[i,2])
        eigvec_range = eigvec_ver[i,:,:]
        eigvec_range = eigvec_range[:,eig_index[i,:]]
        vec1[3*i:3*i+3] = eigvec_range[:,0]
        vec2[3*i:3*i+3] = eigvec_range[:,1]
        vec3[3*i:3*i+3] = eigvec_range[:,2]
        eigvalues[3*i:3*i+3] = eig_ver[i,:]
    newarrayeig.vector().set_local(eigvalues)
    vecarray[0].vector().set_local(vec1)
    vecarray[1].vector().set_local(vec2)
    vecarray[2].vector().set_local(vec3)
    funarray.vector().set_local(phi)
    return 
    