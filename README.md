# multiscale_metamaterial
Multiscale model for Metamaterial calculations

### DESCRIPTION OF THE MODEL
These files are designed to perform metamaterial calculations via a multiscale model, where a large problem is divided into the sum of smaller problems. 
 This means that in order to analyze a complete metamaterial probe, instead of doing the direct FEA calculations with all of its degrees of freedom, where the dimension of the problem might be large enough to overflow the RAM of the computer, the problem is instead divided into smaller problems that reduce the total computational cost. 
 The way in which this is done is by dividing the whole geometry of the probe into smaller regions and these regions are solved independently, therefore, reducing the dimension of the problem an reducing the computational cost.
 In order to apply the boundary conditions to each of the different smaller regions, the original boundary conditions of the whole probe are not enough, so a fast calculation of the displacements of the boundaries of each region is performed. 
 This is achieved by solving the whole probe as if it were an homogeneous solid with averaged properties. This is nothing but a coarse approximation of the real solution that later will be "corrected" thanks to the exact metamaterial calculations of each small region in the solid
 This code enables to know the total energy of the probe as well as the energy in each of the smaller regions. Moreover, the energy and deformation of each bar in the matamaterial is available

In this code, three geometries ("multiscale") are used:

	1- When the probe is solved as an homogeneous solid, it is discretized with "thin hexaedra", which are used for the FEA calculation
 
	2- To analise certain regions of the solid with the exact metamaterial, the solid must be divided into regions, which are also hexaedra but larger than the latter (these are called "large hexaedra") and they are not used for any FEA calculation
 
	3- The different large hexaedra are filled with the corresponding metamaterial bars. This is why they are referred as "metamaterial regions"

The homogeneuos solid calculation is used to impose approximate boundary conditions for each of the metamaterial regions. This is done by relating the interior displacements in the homogeneous solid to the nodes of the metamaterial in the boundaries of each region.
 This coupling ensures that the displacements in adjacent boundaries of two metamaterial regions are the same. However, the sum of the forces in an interior node should be zero (equilibrium in each node) but given that this condition is not imposed, it is not going to be verified.
 To solve this, an optimization method (gradient method) is used to correct this by building an objective function where the sum of the forces in every node is considered. To minimize this function, the position of the nodes is slightly modified. 
 The result of this movement is that the displacements imposed in the boundaries of each metamaterial region are no longer an approximation obtained with the homogeneous solid, but they are the exact displacements of these nodes.
 

### EXECUTION PROCESS
1- Open the file "MAIN_multiescala_hexaedros_mallado_fino_grueso_MEF"

2- Choose the name for the generated results documents

3- Now, three aspects can be modified:

        1) The dimensions of the probe (it is a prism) and its dimension in X, Y and Z can be selected.
	
	2) The boundary conditions (both for the fixed points and for the imposed displacements)
 
               Given that the boundary conditions are imposed to the homogeneous solid, only three degrees of freedom are considered in each node (displacements in X, Y and Z)
	       
               The coordinates where the conditions are imposed may be specific points in space or they may be specified as "all" if all the points with a specific coordinate are to be selected.
               There are different predesigned boundary conditions that are commented, so with the comment/uncomment command one can be selected (tension-compression, bending, torsion)
	       
	3) The hyperparameters to specify how to solve the problem. The hyperparameters range from 7 to 11 depending on the method for solving the problem. 
 
    		a) if the method is gradient: "off", the 7 hyperparameters are:
    			-iter_desplaz = Number of steps in which the total final displacement is divided      
			-iter_max = Number of iterations to converge in each step
			-tolerancia = Convergence critera: difference between the error vector in two succesive iterations (in the same step)
			-tolerancia_2 = Convergence critera: difference between the error vector in two non-succesive iterations (in the same step)
			-diferencial = Longitude differential for the calculation of derivatives
			-eps_limite = Deformation threshold in a hexaedron. If it is passed, the the corresponding region is analyzed via metamaterial
			-modelo_calculo = Model to update the tangent Young's modulus of a given hexaedron (3 possibilities) 
                                                "modelo_energia_hexaedros_gruesos": This model updates all thin hexahedrons with the same value as the thick hexahedron. It updates them by adjusting the tangent Young's modulus so that the energy matches that of the zone calculated by the metamaterial
                                                "modelo_daño_hexaedros_finos": This model updates all thin hexahedrons with different values from the thick hexahedron. It updates them by adjusting the tangent Young's modulus so that the energy of each thin hexaedron matches that calculated by the metamaterial 
                                                "modelo_daño_hexaedros_gruesos": This model updates all thin hexahedrons with the same value as the thick hexahedron. The Young is updated using a damage parameter
			-exponente_reparto_young = Necessary only for "modelo_daño_hexaedros_finos". If the value is 0, then all thin hexaedrons inside the large hexaedron are updated with the same value of tangent Young's modulus
		
		b) if the method is gradient: "on", there are 11 hyperparameters:
			-The first seven ones are the same as before
			-step_gradiente_desplaz = constant for updating displacement values of the boundary nodes of each region
			-step_gradiente_giros = constant for updating rotations in the boundary nodes of each region
			-iter_max_solido_grad = Number of iterations to converge in each step (This value substitutes "iter_max")
			-iter_max_gradiente = Maximum number of iterations to converge using the gradient method
4- Results: 

	a) During the execution of the code
		- Information of the beginning and end of each command is shown
		- Detailed information of the total energy accumulated, increments of energy, deformation, tangent Young's modulus, etc. is shown for each hexaedron
	b) 11 Files generated:
		- "multiescala_energ_element_filename.txt": total energy accumulated (after the different displacement steps) of every large hexaedron calculated with the homogeneous solid model  
		- "multiescala_inc_energ_element_filename.txt": increment energy (in the displacement step) of every large hexaedron calculated with the homogeneous solid model  
		- "multiescala_energ_total_filename.txt": total energy accumulated (after the different displacement steps) in the whole solid with the homogeneous solid model
		- "multiescala_inc_energ_total_filename.txt": increment energy (in the displacement step) in the whole solid with the homogeneous solid model
		- "multiescala_fuerza_filename.txt": total force applied (after the different displacement steps) in the boundary pints where the displacement is imposed
		- "multiescala_energ_zona_metam_filename.txt": total energy accumulated (after the different displacement steps) of every large hexaedron calculated with the metamaterial 
		- "multiescala_inc_energ_zona_metam_filename.txt": increment energy (in the displacement step) of every large hexaedron calculated with the metamaterial 
		- "multiescala_energ_metam_total_filename.txt": total energy accumulated (after the different displacement steps) in the whole solid with the metamaterial 
		- "multiescala_inc_energ_metam_total_filename.txt": increment energy (in the displacement step) in the whole solid with the metamaterial 
		- "multiescala_curva_tension_deformacion_hexaedros_gruesos_filename.txt": The results obtained allow to approximate the stress-strain curve of the complete probe
		- "multiescala_deformacion_barra_max_hexaedros_gruesos_filename.txt": Maximum deformation detected in a bar in each large hexaedron

