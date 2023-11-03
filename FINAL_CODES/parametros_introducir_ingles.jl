

### GENERATE RESULT FILES
numero_fichero="file_2023"   # Name of the file

### PLOT PARAVIEW
plot_paraview=false          # To generate or not the paraview files for visializing the problem


### DIMENSIONS OF THE RECTANGULAR PROBE AND THE METAMATERIAL CELLS
Nx=12   # cell number in x
Ny=12   # cell number in en y
Nz=12   # cell number in en z

X=10*Nx    # x dimension of the probe
Y=10*Ny    # x dimension of the probe
Z=10*Nz    # x dimension of the probe
radio=5    # radius of the metamaterial bars


### NUMBER OF DIVISIONS IN THE DIFFERENT AXES TO GENERATE THE LARGE HEXAEDRA
N_3D_x_grueso=4
N_3D_y_grueso=4
N_3D_z_grueso=4


### NUMBER OF DIVISIONS IN THE DIFFERENT AXES TO GENERATE THE THIN HEXAEDRA (SOLID MESH)
N_3D_x_fino=8
N_3D_y_fino=8
N_3D_z_fino=8


### UNITARY CELL FOR THE METAMATERIAL
tipo_celdilla="celdilla_FCC"    
                                #  "celdilla_simple"
                                #  "celdilla_FCC" 


### HYPERPARAMETERS
gradiente=false         #To apply or not the gradient method to reduce errors when performing the multiscale coupling

iter_desplaz = 10       # Number of steps in which the total final displacement is divided      
iter_max = 2            #Number of iterations to converge in each step
tolerancia = 1e-5       #Convergence critera: difference between the error vector in two succesive iterations (in the same step)
tolerancia_2 = 1e-6     #Convergence critera: difference between the error vector in two non-succesive iterations (in the same step)
diferencial = 1e-8      #Longitude differential for the calculation of derivatives
eps_limite = 1e-12x     #Deformation threshold in a hexaedron. If it is passed, the the corresponding region is analyzed via metamaterial
modelo_calculo = "modelo_energia_hexaedros_gruesos"      #Model to update the tangent Young's modulus of a given hexaedron (3 possibilities) 
                                #    "modelo_energia_hexaedros_gruesos"
                                #    "modelo_daño_hexaedros_finos": 
                                #    "modelo_daño_hexaedros_gruesos": 
exponente_reparto_young = 0 #Necessary only for "modelo_daño_hexaedros_finos"

#if the method is gradient: "on", there are additional hyperparameters:
step_gradiente_desplaz = 0.02*10/(194000*radio^2)              #constant for updating displacement values of the boundary nodes of each region
step_gradiente_giros = 0.02/(10^2)*10^3/(194000*radio^4)       #constant for updating rotations in the boundary nodes of each region
iter_max_solido_grad = 1                                       #Number of iterations to converge in each step (This value substitutes "iter_max")
iter_max_gradiente = 2                                         #Maximum number of iterations to converge using the gradient method


################################### BOUNDARY CONDITIONS #####################################

### RESTRICTIONS 

# restriction 1
punto_restriccion=[0 "all" "all"]      # Coordinates of the node [x,y,z]
gdl_restriccion=[1 1 1]                # degrees of freedom [ux, uy, uz]



### IMPOSED DISPLACEMENTS 

# displacement 1
punto_desplazamiento=[X "all" "all"]    # Coordinates of the node [x,y,z]
valor_desplazamiento=[ 0  0 10e-3 ]     # degrees of freedom [ux, uy, uz]