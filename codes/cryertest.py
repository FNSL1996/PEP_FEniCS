#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 10:44:36 2022

@author: felipe
"""

'''
Poro-elasto-plastic u-p formulation: With TangentOperator
'''

'''
Libraries used
'''

import numpy as np
import dolfin as df
import time
import auxiliary_updated as au
from datetime import datetime

beginning = datetime.now()
start_time = beginning.strftime('%H:%M:%S')

'''
Compiler options
'''
df.parameters["form_compiler"]["cpp_optimize"] = True
df.parameters["form_compiler"]["representation"] = 'uflacs'

'''
Import mesh files
- It's possible to call a file within a different folder, but I need to modify this next few lines...
'''


folder_name = 'cryertest/'
file_name = 'octant'
mesh = df.Mesh(file_name + '.xml')
facets = df.MeshFunction('size_t',mesh,file_name+'_facet_region.xml')
volumes = df.MeshFunction('size_t',mesh,file_name+'_physical_region.xml')

'''
Gravity setting: 1 for gravity. 0 for non gravity
'''

gravity = 9.81 #m/s²
weight = df.Constant(0.0)

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
Hydraulic and mechanical parameters: Haagenson et al (2020)
'''

k = 1.0e-18 # Permeability (m²)
mu = 8.9e-4 # Dynamic viscosity (Pa*s)
phi = 0.05 # Porosity (-)
rho_f = 1000 # Fluid density (kg/m³)
rho_r = 2800 # Rock or solid density (kg/m³)
K = 1.0e10 # Bulk modulus (Pa) E = 1.5e10
G = 6.0e9 # Shear modulus (Pa)
E = 1.5e10 # Young's modulus (Pa)
Kf = (4.4e-10)**-1 # Fluid bulk modulus (Pa)
Ks = (1.0e-11)**-1 # Solid bulk modulus (Pa)
nu = 0.25 # Poisson's ratio (-)
alpha = 1.0 - K/Ks # Biot's coefficient (-) 
M = ((phi/Kf) + (alpha - phi)/Ks)**-1 # Biot's modulus (Pa)

# Elastoplastic parameters (defined to be a poroelastic test, and compare it with the analytical solution)
sig0 = 100e9 # Yield strength
Et = E/2. # Tangent modulus
H = E*Et/(E - Et) # Hardening modulus

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
beta = df.Function(W0, name='Reduccion de stress deviatorico')
N = df.Function(W_stress, name='Normal vector to yield surface')
plas_strain = df.Function(W0, name="Cumulative plastic strain")
q_VM = df.Function(W0, name = 'von Mises Equivalent Stress')
p_stress = df.Function(W0, name = 'Mean Stress')
evol_tr = df.Function(W0, name = 'Volumetric strain')


'''
Output files for facets and volumes
'''

file_sd = df.File(folder_name+"facets.pvd")
file_sd << facets

file_sd = df.File(folder_name+"subdomains.pvd")
file_sd << volumes

'''
Time configuration
- It's possible to use a variable time scheme but modification would be needed in the solver
'''

tf = 70000
nt = 500
dt = tf/nt

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

'''
Boundary conditions:
- This segment is highly modifiable. The following implementation is for Cryer's Problem
- In this formulation, the prescribed flux is set using a Neumann boundary condition
'''

    
### DISPLACEMENT CONDITIONS ###
cb1 = df.DirichletBC(W.sub(1).sub(1), df.Constant(0.0), \
    facets, 1) 
cb2 = df.DirichletBC(W.sub(1).sub(0), df.Constant(0.0), \
    facets, 3) 
cb3 = df.DirichletBC(W.sub(1).sub(2), df.Constant(0.0), \
    facets, 4) 
cb4 = df.DirichletBC(W.sub(0), df.Constant(0.0), facets, 2)
    
BCS = [cb1,cb2,cb3,cb4]

### Prescribed flux ###

pr_flux = df.Constant(0.0)

### Prescribed tension ###

pr_ten = df.Constant(-5e9)

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
    
    momentum_1 = df.inner(au.sigma_elastoplastic(u_tr, p_tr, alpha, K, G, beta),df.grad(u_te))*df.dx 
    momentum_2 = -pr_ten*df.dot(u_te,norm)*df.ds(subdomain_data = facets, subdomain_id = 2)
    
    
    Momentum = momentum_1 + momentum_2
    
    
    return Mass + Momentum

def TangentOperator(p_trial, u_trial, p_te, u_te, dt, beta, N):
    
    
    '''
    Conservation of mass
    '''
    
    MassDG_u = alpha*df.div(u_trial)*p_te*df.dx
    MassDG_p = (1/M)*p_trial*p_te*df.dx + (k/mu)*dt*df.inner(df.grad(p_trial),df.grad(p_te))*df.dx

    '''
    Conservation of momentum
    '''
    
    MomentumDG_u = df.inner(au.sigma_tang(K,G,H,u_trial,beta,N),df.grad(u_te))*df.dx
    MomentumDG_p = -alpha*p_trial*df.inner(df.Identity(3),df.grad(u_te))*df.dx
    
    return MassDG_u + MassDG_p + MomentumDG_u + MomentumDG_p
    
    
    
t = 0.0
con = 0
Nitermax, tol = 200, 1e-8

for n in range(0,nt):
    t += dt
    now = datetime.now()
    start_time = now.strftime("%H:%M:%S")
    print('Iteration Start Time = ',start_time)
    Residual = WeakForm(X_func, p_te,  u_te, p_ini, u_ini,dt)
    Tangent = TangentOperator(p_trial, u_trial, p_te, u_te, dt, beta, N)
    A, Res = df.assemble_system(Tangent,-Residual,BCS)
    nRes0 = Res.norm("l2")
    nRes = nRes0
    Du = df.interpolate(X_i,W)
    Du0 = Du.vector()
    niter = 0
    while nRes/nRes0  > tol and niter < Nitermax:
        df.solve(A, du.vector(), Res, "mumps")
        Du.assign(Du+du)
        X_func.assign(Du) 
        _, USu = df.split(X_func)
        sig, N, beta, dp_, q, _, p_now, evol_now, edev_now = au.sig_correction(USu,sig_old,plas_strain,K,G,sig0,H)
      #  Activate for after displacement BC
      #  for bc in BCS:  
      #      bc.homogenize()
        A, Res = df.assemble_system(Tangent, -Residual, BCS)
        nRes = Res.norm('l2')
        print('Residual: ' + str(nRes))
        niter += 1
        if niter == Nitermax:
            report.write('Max number of iterations have been surpassed!')
            report.write('\n')
        	
    plas_strain.assign(plas_strain+au.local_project(dp_, W0, df.dx))
    
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
    
    p_avg.assign(df.project(plas_strain, P0,solver_type = 'mumps'))
    
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
    
    
    # Update solution at last time step
    X_n.assign(X_func)
    p_ini, u_ini = X_n.split(True)
    
    X_func = df.Function(W)
    
    now = datetime.now()
    start_time = now.strftime("%H:%M:%S")
    print('Iteration End Time = ' + start_time)
    con +=1

    
now = datetime.now()
start_time = now.strftime("%H:%M:%S")
print('End Time = ' + start_time)

