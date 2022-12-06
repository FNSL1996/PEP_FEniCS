#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 10:44:36 2022

@author: felipe
"""

'''
Poro-elasto-plastic u-p formulation
Full Model Run
'''

'''
Libraries used
'''

import numpy as np
import dolfin as df
import time
import auxiliary_updated as au
from datetime import datetime

folder_name = 'SR_1/'
file_name = '../meshes/newmesh'

report = open(folder_name+'Code_Report.txt','w',encoding = 'utf-8')
report.write('Code Report: Geothermal Model Run - Poro-elasto-plastic')
report.write('\n')
report.write('\n')

beginning = datetime.now()
start_time = beginning.strftime('%H:%M:%S')
report.write('Start Time = ' + start_time)
report.write('\n')
report.write('\n')

'''
Compiler options
'''
df.parameters["form_compiler"]["cpp_optimize"] = True
df.parameters["form_compiler"]["representation"] = 'uflacs'

'''
Import mesh files
'''


mesh = df.Mesh(file_name + '.xml')
facets = df.MeshFunction('size_t',mesh,file_name+'_facet_region.xml')
volumes = df.MeshFunction('size_t',mesh,file_name+'_physical_region.xml')

'''
Finite elements configuration
'''

ele_p  = df.FiniteElement("CG",  mesh.ufl_cell(), 1) # Pressure
ele_us  = df.VectorElement("CG",  mesh.ufl_cell(), 2) # Solid displacement
W = df.MixedElement([ele_p, ele_us])
W = df.FunctionSpace(mesh, W)
We = df.TensorElement('DG',mesh.ufl_cell(),1,shape=(3,3)) # Tensor element for stress tensor
W_stress = df.FunctionSpace(mesh, We) # Solution space for stress tensor
W0e = df.FiniteElement("DG", mesh.ufl_cell(), 1) # Finite element for cumulative plastic strain
W0 = df.FunctionSpace(mesh, W0e) # Solution space for cumulative plastic strain

'''
Solution spaces
'''

Z = df.FunctionSpace(mesh, "CG", 2)
V = df.VectorFunctionSpace(mesh, "CG", 1)
S = df.TensorFunctionSpace(mesh, "CG", 1)
P = df.FunctionSpace(mesh, "CG", 1)


'''
Hydraulic and mechanical parameters: Berea Sandstone
'''

k = 1.9e-13 # Permeability (m²)
mu = 1e-3 # Dynamic viscosity (Pa*s)
phi = 0.19 # Porosity (-)
rho_f = 1000 # Fluid density (kg/m³)
rho_r = 2800 # Rock or solid density (kg/m³)
K = 8.0e9 # Bulk modulus (Pa)
G = 6.0e9 # Shear modulus (Pa)
E1 = 9*K*G/(3*K + G) # Young's Modulus (Pa)
Kf = (2.25e9) # Fluid bulk modulus (Pa)
Ks = (3.6e10) # Solid bulk modulus (Pa)
nu = 0.2 # Poisson's ratio (-)
alpha = 1.0 - K/Ks # Biot's coefficient (-) 
M = ((phi/Kf) + (alpha - phi)/Ks)**-1 # Biot's modulus (Pa)
sig0_1 = df.Constant(5.0e6) # Yield strength
Et1 = E1/5. # Tangent modulus
H1 = E1*Et1/(E1 - Et1) # Hardening modulus

'''
Fault' Shear Modulus
'''

G_falla = 6.0e9 # Shear modulus (Pa)´
E2 = 9*K*G_falla/(3*K + G_falla) # Young's Modulus (Pa)
sig0_2 = df.Constant(5.0e6) # Yield strength
Et2 = E2/2.5 # Tangent modulus
H2 = E2*Et2/(E2 - Et2) # Hardening modulus

'''
Trial and Test Functions
'''

X_func = df.Function(W)

p_tr, u_tr = df.split(X_func)

p_trial, u_trial = df.TrialFunctions(W)

p_te, u_te = df.TestFunctions(W)

Du = df.Function(W)
du = df.Function(W)
U = df.Function(W)

'''
Functions ot keep track of the current internal state and increments
'''

sig = df.Function(W_stress)
edev_tensor = df.Function(W_stress, name = 'Deviatoric Strain tensor')
sig_old = df.Function(W_stress)
beta = df.Function(W0, name='Deviatoric stress reduction') 
N = df.Function(W_stress, name='Normal vector to yield surface')
plas_strain = df.Function(W0, name="Cumulative plastic strain")
q_VM = df.Function(W0, name = 'von Mises Equivalent Stress')
p_stress = df.Function(W0, name = 'Mean Stress')
evol_tr = df.Function(W0, name = 'Volumetric strain')

# Domain's functions
sig1 = df.Function(W_stress)
edev_tensor1 = df.Function(W_stress)
sig_old1 = df.Function(W_stress)
beta1 = df.Function(W0, name='Deviatoric stress reduction') 
N1 = df.Function(W_stress, name='Normal vector to yield surface')
plas_strain1 = df.Function(W0, name="Cumulative plastic strain")
q_VM1 = df.Function(W0, name = 'von Mises Equivalent Stress')
p_stress1 = df.Function(W0, name = 'Mean Stress')
evol_tr1 = df.Function(W0, name = 'Volumetric strain')

# Fault's functions
sig2 = df.Function(W_stress)
edev_tensor2 = df.Function(W_stress)
sig_old2 = df.Function(W_stress)
beta2 = df.Function(W0, name='Deviatoric stress reduction') 
N2 = df.Function(W_stress, name='Normal vector to yield surface')
plas_strain2 = df.Function(W0, name="Cumulative plastic strain")
q_VM2 = df.Function(W0, name = 'von Mises Equivalent Stress')
p_stress2 = df.Function(W0, name = 'Mean Stress')
evol_tr2 = df.Function(W0, name = 'Volumetric strain')

'''
Output files for facets and volumes
'''

file_sd = df.File(folder_name+"facets.pvd")
file_sd << facets

file_sd = df.File(folder_name+"subdomains.pvd")
file_sd << volumes

'''
Time configuration
'''

tf = 2.6e6 + 100000*200

dt = 0

'''
Initial conditions
'''

X_i = df.Expression(
        (
            '0.0',       		# p    
            '0.0','0.0','0.0'  	# (us1, us2, us3)
        ),degree=2)
X_n = df.interpolate(X_i,W) # Initial Condition

p_ini, u_ini = X_n.split(True) # Initial Condition for every field

'''
Boundary conditions:
- This segment is highly modifiable. The following implementation is for the Geothermal Synthetic Model
- In this formulation, the prescribed flux is set using a Neumann boundary condition
'''

def BoundaryConditions(t,dp):
    
    ### DISPLACEMENT CONDITIONS ###
    # Bottom
    cb1 = df.DirichletBC(W.sub(1).sub(2), df.Constant(0.0), \
        facets, 16) 
    # South
    cb2 = df.DirichletBC(W.sub(1).sub(1), df.Constant(0.0), \
        facets, 17) 
    # North
    cb3 = df.DirichletBC(W.sub(1).sub(1), df.Constant(0.0), \
        facets, 19)
    # West
    cb4 = df.DirichletBC(W.sub(1).sub(0), df.Constant(0.0), \
        facets, 20)
    # East
    cb5 = df.DirichletBC(W.sub(1).sub(0), df.Constant(0.0), \
        facets, 18)
    # Fault displacement condition
    cb6 = df.DirichletBC(W.sub(1).sub(1), df.Expression('beta',\
                              beta=beta,degree=2,domain=mesh), facets,100)
    cb7 = df.DirichletBC(W.sub(1).sub(1), df.Constant(0.0),\
        facets, 14)
    # Free surface condition
    cb8 = df.DirichletBC(W.sub(0), df.Constant(0.0), facets, 15)
    # Geothermal reservoir overpressure
    cb9 = df.DirichletBC(W.sub(0), df.Constant(dp), facets, 1)
    cb10 = df.DirichletBC(W.sub(0), df.Constant(dp), facets, 2)
    cb11 = df.DirichletBC(W.sub(0), df.Constant(dp), facets, 3)
    cb12 = df.DirichletBC(W.sub(0), df.Constant(dp), facets, 4)
    cb13 = df.DirichletBC(W.sub(0), df.Constant(dp), facets, 5)
    cb14 = df.DirichletBC(W.sub(0), df.Constant(dp), facets, 6)
    cb15 = df.DirichletBC(W.sub(0), df.Constant(dp), facets, 7)
    cb16 = df.DirichletBC(W.sub(0), df.Constant(dp), facets, 8)
    
    if t > 4 + 100000*200:
        return [cb1,cb2,cb3,cb4,cb5,cb6,cb7,cb8,cb9,cb10,cb11,cb12,cb13,cb14,cb15,cb16]
    else: # Pre-slip condition
        return [cb1,cb2,cb3,cb4,cb5,cb8,cb9,cb10,cb11,cb12,cb13,cb14,cb15,cb16]


### Prescribed flux ###

pr_flux = df.Constant(0.0)

### Prescribed tension ###

pr_ten = df.Constant(0.0)

'''
Normal vector and Trial, test, and solution functions 
'''

n = df.FacetNormal(mesh)
norm = df.as_vector([n[0],n[1],n[2]])
X = df.Function(W)

'''
Fluid source/sink
'''

ff = df.Constant(0.0)

'''
Configuration for plastic strain 
'''

P0 = df.FunctionSpace(mesh, "DG", 1)
p_avg = df.Function(P0, name="Accumulated Plastic strain")
p_avg_1 = df.Function(P0, name="Plastic strain")
p_avg_2 = df.Function(P0, name="Plastic strain")
S1 = df.TensorFunctionSpace(mesh, 'DG', 1, shape=(3,3))
W_eigen = df.VectorFunctionSpace(mesh, 'DG', 1)

'''
Weak Form Assembly
'''

def WeakForm(X_func, p_te, u_te, p_ini, u_ini, dt):
        
    p_tr, u_tr = df.split(X_func)
    
    '''
    Conservation of mass
    '''
    
    mass_1 = alpha*df.div(u_tr)*p_te*df.dx
    mass_2 = -alpha*df.div(u_ini)*p_te*df.dx
    mass_3 = (1/M)*p_tr*p_te*df.dx
    mass_4 = -(1/M)*p_ini*p_te*df.dx
    mass_5 = (k/mu)*dt*df.inner(df.grad(p_tr),df.grad(p_te))*df.dx

    Mass = mass_1 + mass_2 + mass_3 + mass_4 + mass_5 
    

    '''
    Conservation of momentum
    '''
    
    momentum_1 = df.inner(au.sigma_elastoplastic(u_tr, p_tr, alpha, K, G, beta1),df.grad(u_te))*df.dx(subdomain_data = volumes, subdomain_id = 3) + \
          df.inner(au.sigma_elastoplastic(u_tr, p_tr, alpha, K, G_falla, beta2),df.grad(u_te))*df.dx(subdomain_data = volumes, subdomain_id = 2)
    
    Momentum = momentum_1
    
    return Mass + Momentum

def TangentOperator(p_trial, u_trial, p_te, u_te, dt, H1, beta1, N1, H2, beta2, N2):
    
    
    '''
    Conservation of mass
    '''
    
    MassDG_u = alpha*df.div(u_trial)*p_te*df.dx
    MassDG_p = (1/M)*p_trial*p_te*df.dx + (k/mu)*dt*df.inner(df.grad(p_trial),df.grad(p_te))*df.dx

    '''
    Conservation of momentum
    '''
    
    MomentumDG_u = df.inner(au.sigma_tang(K,G,H1,u_trial,beta1,N1),df.grad(u_te))*df.dx(subdomain_data = volumes, subdomain_id = 3) + \
                  df.inner(au.sigma_tang(K,G_falla,H2,u_trial,beta2,N2),df.grad(u_te))*df.dx(subdomain_data = volumes, subdomain_id = 2)   
    MomentumDG_p = -alpha*p_trial*df.inner(df.Identity(3),df.grad(u_te))*df.dx
    
    return MassDG_u + MassDG_p + MomentumDG_u + MomentumDG_p
    
    
    
t = 0.0
con = 0
L = 10000
H = 5000
beta = 0
slip_rate = 1.0
dp = 0
n = 0
Nitermax, tol = 400, 1e-8
iters = []
while t < tf:
    if t < 100000*100:
        dt = 100000*20
        dp += 0.05*20*1e6 
    elif t >= 100000*100 and t < 100000*200:
        dt = 200000*10       
    elif t >= 100000*200 and t < 4 + 100000*200:
        dt = 0.5
    elif t >= 4 + 100000*200 and t < 6 + 100000*200:
        dt = 0.25
        beta += dt*slip_rate
    elif t >= 6 + 100000*200 and t < 60 + 100000*200:
        dt = 0.5
    elif t >= 60 + 100000*200 and t < 3600 + 100000*200:
        dt = 60
    elif t >= 3600 + 100000*200 and t < 3600*24 + 100000*200:
        dt = 1800
    else:
        dt = 3600*24 
    t += dt
    now = datetime.now()
    start_time = now.strftime("%H:%M:%S")
    report.write('Iteration Start Time = ' + start_time)
    report.write('\n')
    Residual = WeakForm(X_func, p_te,  u_te, p_ini, u_ini,dt)
    Tangent = TangentOperator(p_trial, u_trial, p_te, u_te, dt, H1, beta1, N1, H2, beta2, N2)
    BCS = BoundaryConditions(t,dp)
    A, Res = df.assemble_system(Tangent,-Residual,BCS)
    nRes0 = Res.norm("l2")
    nRes = nRes0
    Du = df.interpolate(X_i,W)
    Du0 = Du.vector()
    report.write('Increment: ' + str(n + 1))
    report.write('\n')
    niter = 0
    while nRes/nRes0  > tol and niter < Nitermax:
        df.solve(A, du.vector(), Res, "mumps")
        Du.assign(Du+du)
        X_func.assign(Du) 
        _, USu = df.split(X_func)
        sig_1, N_1, beta_1, dp_1, q_1, _, p_now_1, evol_now1, edev_now1 = au.sig_correction(USu,sig_old1,plas_strain1,K,G,sig0_1,H1)
        sig_2, N_2, beta_2, dp_2, q_2, _, p_now_2, evol_now2, edev_now2 = au.sig_correction(USu,sig_old2,plas_strain2,K,G_falla,sig0_2,H2)
        # Local projection for host rock
        au.local_project(sig_1, W_stress, df.dx(subdomain_data = volumes, subdomain_id = 3),sig1)
        au.local_project(edev_now1, W_stress, df.dx(subdomain_data = volumes, subdomain_id = 3),edev_tensor1)
        au.local_project(N_1, W_stress, df.dx(subdomain_data = volumes, subdomain_id = 3), N1)
        au.local_project(beta_1, W0, df.dx(subdomain_data = volumes, subdomain_id = 3), beta1)
        au.local_project(q_1, W0, df.dx(subdomain_data = volumes, subdomain_id = 3), q_VM1)
        au.local_project(p_now_1, W0, df.dx(subdomain_data = volumes, subdomain_id = 3), p_stress1)
        au.local_project(evol_now1, W0, df.dx(subdomain_data = volumes, subdomain_id = 3), evol_tr1)
        # Local projection for fault
        au.local_project(sig_2, W_stress, df.dx(subdomain_data = volumes, subdomain_id = 2),sig2)
        au.local_project(edev_now2, W_stress, df.dx(subdomain_data = volumes, subdomain_id = 2),edev_tensor2)
        au.local_project(N_2, W_stress, df.dx(subdomain_data = volumes, subdomain_id = 2), N2)
        au.local_project(beta_2, W0, df.dx(subdomain_data = volumes, subdomain_id = 2), beta2)
        au.local_project(q_2, W0, df.dx(subdomain_data = volumes, subdomain_id = 2), q_VM2)
        au.local_project(p_now_2, W0, df.dx(subdomain_data = volumes, subdomain_id = 2), p_stress2)
        au.local_project(evol_now2, W0, df.dx(subdomain_data = volumes, subdomain_id = 2), evol_tr2)
        for bc in BCS:  
            bc.homogenize()
        A, Res = df.assemble_system(Tangent, -Residual, BCS)
        nRes = Res.norm('l2')
        report.write('Residual: ' + str(nRes))
        report.write('\n')
        niter += 1
        if niter == Nitermax:
            report.write('Max number of iterations have been surpassed!')
            report.write('\n')

    # Combine vectors for paraview and postprocessing
    # Note that, depending on the specific version of FEniCS you are using, the way in which vectors are assemble may change
    # Because of this, directly adding the vector may not work. combine_results in auxiliary_updated module may solve this problem
    sig.vector().set_local(sig1.vector()[:]+sig2.vector()[:])
    edev_tensor.vector().set_local(edev_tensor1.vector()[:]+edev_tensor2.vector()[:])
    evol_tr.vector().set_local(evol_tr1.vector()[:]+evol_tr2.vector()[:])
    q_VM.vector().set_local(q_VM1.vector()[:]+q_VM2.vector()[:])
    p_stress.vector().set_local(p_stress1.vector()[:]+p_stress2.vector()[:])
        	
    plas_strain1.assign(plas_strain1+au.local_project(dp_1, W0, df.dx(subdomain_data = volumes, subdomain_id = 3)))
    plas_strain2.assign(plas_strain2+au.local_project(dp_2, W0, df.dx(subdomain_data = volumes, subdomain_id = 2)))
    
    iters.append(niter)
    np.savetxt(folder_name+'Iterations_'+str(con)+'.txt',np.array([niter]))
    
    # Postprocessing
    
    p, us = X_func.split(True)
    
    p.rename('Fluid Pressure [Pa]','Presion de poros')
    us.rename('Solid displacement [m]','Desplazamiento del esqueleto solido')
    
    file = df.File(folder_name+"FluidPressure_"+str(con)+".pvd")
    file << p
    file = df.File(folder_name+"US_"+str(con)+".pvd")
    file << us   
    
    flux = -(k/mu)*df.grad(p)
    flux = df.project(flux,V,solver_type = 'mumps')
    flux.rename('Darcy velocity [m/s]','Velocidad de Darcy')
    file = df.File(folder_name+"Flux"+str(con)+".pvd")
    file << flux 
    
    p_avg_1.assign(df.project(plas_strain1, P0,solver_type = 'mumps'))
    p_avg_2.assign(df.project(plas_strain2, P0,solver_type = 'mumps'))
    
    p_avg.vector().set_local(au.combine_results(p_avg_1.vector()[:],p_avg_2.vector()[:]))
    
    file = df.File(folder_name+"pstrain_test"+str(con)+".pvd")
    file << p_avg    
       
    sigma = df.project(sig,S1,solver_type = 'mumps')   
    sigma.rename('Cauchy stress tensor [Pa]','Tensor de Cauchy total') 
    
    vonMises = df.project(q_VM,P0,solver_type = 'mumps')
    vonMises.rename('von Mises Stress [Pa]','Esfuerzo de von Mises') 
    
    meanstress = df.project(p_stress,P0,solver_type = 'mumps')
    meanstress.rename('Mean Stress [Pa]','Esfuerzo medio') 
    
    dev_strain = df.project(edev_tensor, S1,solver_type = 'mumps')
    dev_strain.rename('Deviatoric Strain Tensor [-]', 'Tensor de strain deviatorico')
    
    vol_strain = df.project(evol_tr, P0,solver_type = 'mumps')
    vol_strain.rename('Volumetric Strain [-]', 'Strain volumetrico')
    
    file = df.File(folder_name+"Sigma_"+str(con)+".pvd")
    file << sigma
    file = df.File(folder_name+"VM_"+str(con)+".pvd")
    file << vonMises
    file = df.File(folder_name+"MeanStress_"+str(con)+".pvd")
    file << meanstress  
    file = df.File(folder_name+"DevStrain_"+str(con)+".pvd")
    file << dev_strain 
    file = df.File(folder_name+"VolStrain_"+str(con)+".pvd")
    file << vol_strain      
    
    # Principal stresses, directions, and phi
    
    sigma_eig = df.Function(W_eigen)
    sigma_eig.rename('Principal Stresses [MPa]','Principal Stresses [MPa]')
    eig_1 = df.Function(W_eigen)
    eig_1.rename('Principal Direction 1 [-]', 'Principal Direction 1 [-]')
    eig_2 = df.Function(W_eigen)
    eig_2.rename('Principal Direction 2 [-]', 'Principal Direction 2 [-]')
    eig_3 = df.Function(W_eigen)
    eig_3.rename('Principal Direction 3 [-]', 'Principal Direction 3 [-]')
    phi_tensor = df.Function(W0)
    phi_tensor.rename('Stress tensor shape ratio [-]', 'Stress tensor shape ratio [-]')
    
    au.EigenSolver_DG_Order1(mesh,sig,sigma_eig,[eig_1,eig_2,eig_3],phi_tensor)
    
    file = df.File(folder_name+"EigSigmaTotal_"+str(con)+".pvd")
    file << sigma_eig
    
    file = df.File(folder_name+"Dir1Total_"+str(con)+".pvd")
    file << eig_1
    
    file = df.File(folder_name+"Dir2Total_"+str(con)+".pvd")
    file << eig_2
    
    file = df.File(folder_name+"Dir3Total_"+str(con)+".pvd")
    file << eig_3
    
    file = df.File(folder_name+"PhiTotal_"+str(con)+".pvd")
    file << phi_tensor
    
    # Update solution at last time step
    X_n.assign(X_func)
    p_ini, u_ini = X_n.split(True)
    
    X_func = df.Function(W)
    
    now = datetime.now()
    start_time = now.strftime("%H:%M:%S")
    report.write('Iteration End Time = ' + start_time)
    report.write('\n')
    report.write('\n')
    con +=1
    n += 1

    
now = datetime.now()
start_time = now.strftime("%H:%M:%S")
report.write('End Time = ' + start_time)
report.close()
np.savetxt(folder_name+'NumberIterations.txt',np.array(iters))
