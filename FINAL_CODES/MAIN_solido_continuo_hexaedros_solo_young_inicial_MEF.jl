include(raw"representar_paraview.jl"); using .representar_paraview
include(raw"representar_paraview_sucesion.jl"); using .representar_paraview_sucesion
include(raw"mallado_solido_tridimensional_hexaedros.jl"); using .mallado_solido_tridimensional_hexaedros
include(raw"curva_sigma_epsilon_conjunto.jl"); using .curva_sigma_epsilon_conjunto
include(raw"funciones_MEF_solido_hexaedros.jl"); using .funciones_MEF_solido_hexaedros

using Plots


### DIMENSIONES PROBETA RECTANGULAR Y CELDILLAS METAMATERIAL 
p=20
X=10*p   #dimension x probeta
Y=10*p #dimension y probeta
Z=10*p #dimension z probeta


### NUMERO DE DIVISIONES EN LOS DIFERENTES EJES PARA GENERAR EL MALLADO DEL SOLIDO TRIDIMENSIONAL
q=10
N_3D_x=q
N_3D_y=q
N_3D_z=q


coords_3D, conectividad_3D=mallado_3D_hexaedros(X,Y,Z,N_3D_x,N_3D_y,N_3D_z)
volumenes_hexaedros= volumenes_elementos_hexaedros(coords_3D, conectividad_3D)
println("geometria generada")


### HIPERPARAMETROS
iter_max=3   #Numero de iteraciones hasta alcanzar el valor del desplazamiento final


### CONDICIONES DE CONTORNO DE LA PROBETA (Puntos de la geometria donde restringimos o imponemos un desplazamiento): 
    # Para cada restriccion hay que indicar el punto donde se aplica y los GDL restringidos
        # coordenada --> restringir el movimiento de los puntos con esa coordenada en el eje correspondiente; 
        # "all" --> restringir el movimiento a todos los puntos con esas coordenadas
        # GDL restringidos en coordenadas globales: (ux,uy,uz,thetax,thetay,thetaz): 0 implica que no se restringe y 1 que si se restringe
    # Las nuevas restricciones se añaden con su nombre con [punto_restriccion;punto]+[gdl_restriccion;gdl]
    # (Igual para los desplazamientos)


########################### DISTINTAS RESTRICCIONES

### RESTRICCIONES PARA TRACCION, TORSION, FLEXION PURA
# restriccion 1
punto_restriccion=[0 "all" "all"]#[0 "all" "all"] # Coordenadas del nodo [x,y,z]
gdl_restriccion=[1 1 1]#[1 0 0]# # # grados de libertad: [ux, uy, uz]


### RESTRICCIONES PARA TRACCION SIN EFECTOS DE LA BASE
# # restriccion 1
# punto_restriccion=[0 "all" "all"]
# gdl_restriccion=[1 0 0]
# # restriccion 2
# punto=[0 Y/2 "all"]
# gdl=[0 1 0]#[0 0 1]
# punto_restriccion=vcat(punto_restriccion,punto)
# gdl_restriccion=vcat(gdl_restriccion,gdl)
# # restriccion 3
# punto=[0 "all" Z/2]
# gdl=[0 0 1]#[0 1 0]
# punto_restriccion=vcat(punto_restriccion,punto)
# gdl_restriccion=vcat(gdl_restriccion,gdl)


### RESTRICCIONES PARA FLEXION A 3 PUNTOS
# # restriccion 1
# punto_restriccion=[0 "all" "all"]
# gdl_restriccion=[1 1 1]
# # restriccion 2
# punto= [X "all" "all"]
# gdl=[1 1 1]
# punto_restriccion=vcat(punto_restriccion,punto)
# gdl_restriccion=vcat(gdl_restriccion,gdl)


########################### DISTINTOS TIPOS DE CARGA

### TRACCION 
# # desplazamiento 1
# punto_desplazamiento=[X "all" "all"]#
# valor_desplazamiento=[1.5e-6 0  0]


### CARGADO EN UN EXTREMO Y EMPOTRADO EN OTRO
# desplazamiento 1
punto_desplazamiento=[X "all" "all"]#
valor_desplazamiento=[ 0  0 0.9e-0]


### FLEXION PURA (UN GIRO)
# # desplazamiento 1
# punto_desplazamiento=[X "all" Z]# [X "all" "all"]# 
# valor_desplazamiento=[1.5e-6 0  0 ]#[ 0  0 1.5e-6]#
# # desplazamiento 2
# punto=[X "all" 0]
# valor=[-1.5e-6 0 0 ]
# punto_desplazamiento=vcat(punto_desplazamiento,punto)
# valor_desplazamiento=vcat(valor_desplazamiento,valor)

# ## Para el ensayo de flexion pura, siempre que N_3D_z sea par
# for i in 1:N_3D_z/2-1
#     global punto, valor, punto_desplazamiento, valor_desplazamiento
#     punto=[X "all" Z/N_3D_z*i+Z/2]
#     valor=[1.5e-6*2*i/N_3D_z  0  0 ]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
#     punto=[X "all" -Z/N_3D_z*i+Z/2]
#     valor=[-1.5e-6*2*i/N_3D_z  0  0]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
# end


### TORSION
# # desplazamiento 1
# punto_desplazamiento=[X "all" Z]#[X "all" "all"]#
# valor_desplazamiento=[0 -1.5e-6  0 ]#[0 0 0  0 1e-6 0]#
# # desplazamiento 2
# punto=[X "all" 0]
# valor=[0 1.5e-6  0]#
# punto_desplazamiento=vcat(punto_desplazamiento,punto)
# valor_desplazamiento=vcat(valor_desplazamiento,valor)
# # desplazamiento 3
# punto=[X 0 "all"]
# valor=[ 0 0 -1.5e-6]#
# punto_desplazamiento=vcat(punto_desplazamiento,punto)
# valor_desplazamiento=vcat(valor_desplazamiento,valor)
# # desplazamiento 4
# punto=[X Y "all"]
# valor=[ 0 0 1.5e-6]#
# punto_desplazamiento=vcat(punto_desplazamiento,punto)
# valor_desplazamiento=vcat(valor_desplazamiento,valor)

# # para el ensayo de torsion, siempre que N_3D_z y N_3D_y sean par
# for i in 1:N_3D_z/2-1
#     global punto, valor, punto_desplazamiento, valor_desplazamiento
#     punto=[X "all" Z/N_3D_z*i+Z/2]
#     valor=[ 0 -1.5e-6*2*i/N_3D_z 0 ]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
#     punto=[X "all" -Z/N_3D_z*i+Z/2]
#     valor=[ 0 1.5e-6*2*i/N_3D_z 0]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
# end
# for i in 1:N_3D_y/2-1
#     global punto, valor, punto_desplazamiento, valor_desplazamiento
#     punto=[X Y/N_3D_y*i+Y/2 "all" ]
#     valor=[ 0 0 1.5e-6*2*i/N_3D_y ]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
#     punto=[X -Y/N_3D_y*i+Y/2 "all"]
#     valor=[0  0 -1.5e-6*2*i/N_3D_y]#
#     punto_desplazamiento=vcat(punto_desplazamiento,punto)
#     valor_desplazamiento=vcat(valor_desplazamiento,valor)
# end


### FLEXION A 3 PUNTOS (con N_3D_x par)
# # desplazamiento 1
# punto_desplazamiento=[X/2 "all" "all"]
# valor_desplazamiento=[0 0  -1.5e-6 ]



# para llegar al desplazamiento final de forma incremental
valor_desplazamiento=valor_desplazamiento/iter_max



### inicializar variables
N_elementos=length(conectividad_3D[:,1]); N_nodos=length(coords_3D[:,1])
epsilon_voigt_acumulado=zeros(6,8,N_elementos); u_solido=zeros(3*N_nodos); f_solido=zeros(3*N_nodos)
Youngs_tangentes_elementos=curva_tension_deformacion_conjunto(10)*ones(N_elementos)
energia_acumulada=0


for iter_desplazamiento in 1:iter_max

    global u_solido, f_solido, f_elementos, energia_acumulada, epsilon_voigt_acumulado, deform_equivalente, inc_deform_equivalente, inc_energias_elementos, energias_hexaedros_acumulado, f_acumulado_hexaedros, Youngs_tangentes_elementos, fuerza_reaccion_x_iters, fuerza_reaccion_y_iters, fuerza_reaccion_z_iters, desplaz_iters

    println(); println("###################### ITERACION ", iter_desplazamiento ," ##########################")    


    ### gdls de las condiciones de contorno
    gdl_restring_globales_solido, valor_restricciones_globales_solido=condic_contorno_solido(punto_restriccion,gdl_restriccion,punto_desplazamiento,valor_desplazamiento,coords_3D)


    ### obtener incrementos de desplazamientos, fuerzas, ...
    inc_u_solido, inc_f_solido, inc_energias_elementos, deform_equivalente, inc_deform_equivalente, epsilon_voigt_acumulado=resolucion_solido_3D(coords_3D, conectividad_3D, Youngs_tangentes_elementos, epsilon_voigt_acumulado, gdl_restring_globales_solido, valor_restricciones_globales_solido)
    # for i in 1:N_elementos
    #     println(i)
    #     println(deform_equivalente[i]," ",inc_deform_equivalente[i]," ",volumenes_tetraedros[i]); println(inc_energias_elementos[i]," ", 1/2*inc_deform_equivalente[i]^2*Youngs_elementos[i]*volumenes_tetraedros[i])
    # end
    # Energia de los elementos del solido y del solido completo

    energia_acumulada=energia_acumulada+f_solido'*inc_u_solido+sum(inc_energias_elementos)
    println("energia= ", energia_acumulada)


    # desplazamientos y fuerzas totales como suma de los incrementos 
    u_solido=u_solido+inc_u_solido
    f_solido=f_solido+inc_f_solido


    # representar paraview
    plot_paraview_solido_sucesion_hexaedros(coords_3D, conectividad_3D, u_solido, f_solido, iter_desplazamiento)



    ### COMPARAR LOS DISTINTOS INVARIANTES
    #display(plot(inv_1_solido,label="inv_1"));readline();display(plot!(sqrt.(abs.(inv_2_solido)),label="inv_2"));readline();display(plot!(cbrt.(inv_3_solido),label="inv_3"));readline();display(plot!(deform_equivalente,label="def_equiv"))


    ### OBTENER LA FUERZA APLICADA EN LOS PUNTOS CON DESPLAZAMIENTO IMPUESTO
    fuerza_reaccion_x=[0.0]; fuerza_reaccion_y=[0.0]; fuerza_reaccion_z=[0.0]

    # indice=findall(x->(x<X+0.0001)&&(x>X-0.001),coords_3D[:,1])
    # nodos=coords_3D[indice,4]
    # for nodo in nodos
    #     push!(fuerza_x,f_solido[3*nodo-2])
    #     push!(fuerza_y,f_solido[3*nodo-1])
    #     push!(fuerza_z,f_solido[3*nodo])
    # end
    indices=findall(!in(0),valor_restricciones_globales_solido)
    for gdl in gdl_restring_globales_solido[indices]
        if gdl%3==1
            push!(fuerza_reaccion_x,f_solido[gdl])
        elseif gdl%3==2
            push!(fuerza_reaccion_y,f_solido[gdl])
        elseif gdl%3==0
            push!(fuerza_reaccion_z,f_solido[gdl])
        else 
        end
    end

    fuerza_x_total=sum(fuerza_reaccion_x);fuerza_y_total=sum(fuerza_reaccion_y);fuerza_z_total=sum(fuerza_reaccion_z)
    println("f_aplicada_x= ", fuerza_x_total,"  f_aplicada_y= ", fuerza_y_total,"  f_aplicada_z= ", fuerza_z_total)


    ### COMPROBAR QUE EL SUMATORIO DE FUERZAS Y MOMENTOS ES IGUAL A CERO
    println("Comprobar sumatorio fuerzas nulo: ",sum(f_solido[1:3:end])," ",sum(f_solido[2:3:end])," ",sum(f_solido[3:3:end]))    
    # for i in 1:length(coords_3D[:,1])
    #     println(coords_3D[i,1:3],u_solido[3*i-2:3*i])
    # end
    #println(u_solido)


end

