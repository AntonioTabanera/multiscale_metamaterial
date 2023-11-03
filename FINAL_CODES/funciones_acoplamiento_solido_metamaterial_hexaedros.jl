#MODULO PARA PASAR LOS VALORES DE LOS DESPLAZAMIENTOS DEL SOLIDO A LOS DEL METAMATERIAL (EL ACOPLAMIENTO DE LAS ESCALAS)
module funciones_acoplamiento_solido_metamaterial_hexaedros

include(raw"celdillas_metamaterial.jl"); using .celdillas_metamaterial 
include(raw"funciones_MEF_solido_hexaedros.jl"); using .funciones_MEF_solido_hexaedros
include(raw"funciones_MEF_metamaterial.jl"); using .funciones_MEF_metamaterial
include(raw"curva_sigma_epsilon_conjunto.jl"); using .curva_sigma_epsilon_conjunto
include(raw"representar_paraview_sucesion.jl"); using .representar_paraview_sucesion
include(raw"mallado_solido_tridimensional_hexaedros.jl"); using .mallado_solido_tridimensional_hexaedros
include(raw"curva_sigma_epsilon_material.jl"); using .curva_sigma_epsilon_material
using LinearAlgebra
using StaticArrays
using Plots

export nodos_metamaterial_en_elementos, nodos_metamaterial_en_elementos_sin_repetir, clasificacion_gdl_zona_relajacion, clasificacion_gdl_multiescala, clasificacion_gdl_zona, resolucion_iterativa_acoplamiento_zonas_mallado_fino_grueso, clasificar_nodos_metam_en_hexaedros, resolucion_iterativa_acoplamiento_zonas_mallado_fino_grueso_gradiente_fuerzas


# funcion para clasificar los nodos de metamaterial dentro de los hexaedros
# En este caso, las barras que se encuentran entre 2 elementos, se consideran en ambos 
function nodos_metamaterial_en_elementos(X,Y,Z, coords_metamaterial, conectividad_metamaterial, coords_solido, conectividad_solido, volumenes_hexaedros, radio)
  
    N_elementos_solidos=length(conectividad_solido[:,1]) # total de elementos hexaedricos
    N_nodos_metamaterial=length(coords_metamaterial[:,1]) # total de nodos en la geometria del metamaterial
    N_elementos_metam=length(conectividad_metamaterial[:,1])
    #indice_coords_pto_corte=N_nodos_metamaterial
    indice_conectiv_pto_corte=N_elementos_metam
    epsilon=1e-10; young_material=Young_material()
    
    conectiv_metam_hexaedros_exacto=Dict(); conectiv_metam_hexaedros_extendido=Dict(); conectiv_metam_hexaedros_reducido=Dict()  # conectividad de elementos que tienen los dos nodos en el interior del hexaedro o en su frontera
    porcentajes_barras_dentro_hexaedros=Dict()   
    coords_metam_hexaedros_exacto=Dict(); coords_metam_hexaedros_extendido=Dict(); coords_metam_hexaedros_reducido=Dict()
    coords_metam_fronteras_exactas_hexaedros=Dict()
    tensiones_pandeo_elementos=Dict()


    for elemento_solido in 1:N_elementos_solidos

        barras_hexaedro_exacto=[];barras_hexaedro_extendido=[];barras_hexaedro_reducido=[]
        porcentajes_barras_dentro_hexaedro=[]
        porcentajes_barras_dentro_hexaedro_dict=Dict()                       
        Nodos_metam_hexaedro_exacto=[];Nodos_metam_hexaedro_extendido=[];Nodos_metam_hexaedro_reducido=[] 
        Nodos_metam_frontera_exacta_hexaedro=[];Nodos_metam_frontera_extendido_hexaedro=[]



        coords_vertices_hexaedro=zeros(8,3) #4 nodos con 3 coordenadas por nodo (vertices del hexaedro)
        for i in 1:8 
            nodo_solido=conectividad_solido[elemento_solido,i] # numero de nodo de los vertices del solido
            coords_vertices_hexaedro[i,:]=coords_solido[nodo_solido,1:3]  # coordenadas de los nodos del elemento hexaedrico
        end
        volumen_hexaedro = volumenes_hexaedros[elemento_solido] #volumen del hexaedro en el que estamos

        clasificacion_caras_hexaedro=tipo_caras_hexaedro(coords_vertices_hexaedro,X,Y,Z)
        #println(clasificacion_caras_hexaedro)


        # ESTUDIAR CADA NODO DE METAMATERIAL
        for nodo_1 in 1:N_nodos_metamaterial
            coords_nodo_1=coords_metamaterial[nodo_1,1:3]
            valor_1, caras_pertenece_nodo_1=is_in_hexaedro(coords_nodo_1, coords_vertices_hexaedro, volumen_hexaedro)

            # SI ESTA EN EL HEXAEDRO (INTERIOR O CARAS)
            if  valor_1 != 0 #para ver si esta dentro o no se usan y comparan los volumenes
                
                # buscar con que nodos tiene conectividad el nodo 1
                indice_1=findall(in(nodo_1),conectividad_metamaterial[:,1:2])
                N_nodos_contacto=length(indice_1) # numero total de conectividades del nodo

                # hay que asegurarse de que al menos alguna conectividad sea con otro nodo de dentro para que no este aislado
                nodos_2=zeros(Int64,N_nodos_contacto); valores_2=zeros(Int64,N_nodos_contacto); coords_nodos_2=zeros(N_nodos_contacto,3); caras_pertenecen_nodos_2=Dict()

                for i in 1:N_nodos_contacto
                    # obtenemos para cada nodo con el que tiene conectividad, cual es ese nodo
                    if indice_1[i][2]==1
                        nodos_2[i]=conectividad_metamaterial[indice_1[i][1],2]
                    else
                        nodos_2[i]=conectividad_metamaterial[indice_1[i][1],1]
                    end
                    coords_nodos_2[i,:]=coords_metamaterial[nodos_2[i],1:3]
                    valores_2[i], caras_pertenecen_nodos_2[i]=is_in_hexaedro(coords_nodos_2[i,:], coords_vertices_hexaedro, volumen_hexaedro)
                end
                

                # SI NO CONECTA CON NINGUN ELEMENTO DE DENTRO: no consideramos el nodo
                if valores_2==zeros(N_nodos_contacto)

                # SI CONECTA CON ALGUNO DE DENTRO
                else
                    push!(Nodos_metam_hexaedro_exacto, nodo_1); push!(Nodos_metam_hexaedro_reducido, nodo_1)

                    for i in 1:N_nodos_contacto

                        valor_2=valores_2[i]
                        nodo_2=nodos_2[i]
                        coords_nodo_2=coords_nodos_2[i,:]
                        caras_pertenece_nodo_2=caras_pertenecen_nodos_2[i]


                        # SI EL NODO 1 ESTA EN EL INTERIOR
                        if caras_pertenece_nodo_1 == [] # en este caso estará en el interior del hexaedro
                            
                            # SI EL NODO 2 TAMBIEN ESTA EN EL INTERIOR: se cuenta el elemento entero
                            if valor_2 == 1              
                                barra=conectividad_metamaterial[indice_1[i][1],3]    
                                push!(barras_hexaedro_exacto,barra); push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                push!(porcentajes_barras_dentro_hexaedro, 1)

                            # SI EL NODO 2 ESTA EN LAS CARAS: se cuenta el elemento entero
                            elseif valor_2 == 2
                                barra=conectividad_metamaterial[indice_1[i][1],3] 
                                push!(barras_hexaedro_exacto, barra); push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                push!(porcentajes_barras_dentro_hexaedro, 1)
                                push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)

                            # SI EL NODO 2 ESTA FUERA: se obtiene el punto de corte y se genera un nuevo elemento
                            else
                                # aumentamos las matrices de coordenadas (si el nuevo nodo no está ya metido) y de conectividad de metamaterial con el nodo del punto de corte
                                coords_punto_corte= punto_corte_barras_caras_hexaedro(coords_nodo_1,coords_nodo_2,coords_vertices_hexaedro,volumen_hexaedro)[1]

                                # hay que ver si el nodo esta ya metido o no. Para ello, buscamos solo entre las coordenadas nuevas
                                k=false; p=0
                                for indice in N_nodos_metamaterial:length(coords_metamaterial[:,1]) 
                                    if abs(coords_punto_corte[1]-coords_metamaterial[indice,1])<epsilon
                                        if abs(coords_punto_corte[2]-coords_metamaterial[indice,2])<epsilon
                                            if abs(coords_punto_corte[3]-coords_metamaterial[indice,3])<epsilon
                                                p=indice
                                                k=true; break   
                                            else
                                            end
                                        else
                                        end
                                    else
                                    end
                                end

                                indice_conectiv_pto_corte=indice_conectiv_pto_corte+1
                                # si no se ha registrado el nodo todavia, le registramos ahora
                                if k==false
                                    indice_coords_pto_corte=length(coords_metamaterial[:,1])+1
                                    coords_metamaterial=vcat(coords_metamaterial,[coords_punto_corte indice_coords_pto_corte])
                                    conectividad_metamaterial=vcat(conectividad_metamaterial,[nodo_1 indice_coords_pto_corte indice_conectiv_pto_corte])
                                # si ya se habia registrado el nodo
                                else
                                    indice_coords_pto_corte=p
                                    conectividad_metamaterial=vcat(conectividad_metamaterial,[nodo_1 indice_coords_pto_corte indice_conectiv_pto_corte])
                                end
                               

                                push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                push!(porcentajes_barras_dentro_hexaedro,1)
                                push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                #push!(Nodos_metam_hexaedro_reducido, nodo_2)
                            end


                        # SI EL NODO 1 ESTA EN ALGUNA CARA
                        else 
                            # SI EL NODO 2 ESTA EN EL INTERIOR: se cuenta el elemento entero
                            if valor_2 == 1
                                barra=conectividad_metamaterial[indice_1[i][1],3] 
                                push!(barras_hexaedro_exacto, barra); push!(barras_hexaedro_extendido, barra); push!(barras_hexaedro_reducido, barra)
                                push!(porcentajes_barras_dentro_hexaedro, 1)
                                push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)  
                                #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                            # SI EL NODO 2 ESTA EN LAS CARAS: se cuenta parte del elemento
                            elseif valor_2 == 2

                                # Hay que determinar cuantas caras tienen en comun y cuales son
                                caras_comunes=[] 
                                for cara in caras_pertenece_nodo_1 
                                    if findall(in(cara), caras_pertenece_nodo_2) != []
                                        push!(caras_comunes, cara)
                                    else
                                    end
                                end
                                N_caras_comunes=length(caras_comunes)

                                # SI NO TIENEN NINGUNA CARA EN COMUN
                                if N_caras_comunes==0
                                    barra=conectividad_metamaterial[indice_1[i][1],3] 
                                    push!(barras_hexaedro_exacto, barra);push!(barras_hexaedro_extendido, barra); push!(barras_hexaedro_reducido, barra)
                                    push!(porcentajes_barras_dentro_hexaedro, 1)
                                    push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                    push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                    #push!(Nodos_metam_hexaedro_reducido, nodo_2)


                                # SI SOLO TIENEN UNA CARA EN COMUN
                                elseif N_caras_comunes==1

                                    # SI LA CARA COMUN ES DEL EXTERIOR DEL SOLIDO: se cuenta entero el elemento
                                    if clasificacion_caras_hexaedro[caras_comunes] == [1]
                                        barra=conectividad_metamaterial[indice_1[i][1],3] 
                                        push!(barras_hexaedro_exacto, barra); push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                        push!(porcentajes_barras_dentro_hexaedro, 0.5) #1
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                    # SI LA CARA COMUN NO ES DEL EXTERIOR DEL SOLIDO: se cuenta la mitad
                                    else
                                        barra=conectividad_metamaterial[indice_1[i][1],3] 
                                        push!(barras_hexaedro_exacto, barra);push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                        push!(porcentajes_barras_dentro_hexaedro, 0.5)
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        #push!(Nodos_metam_hexaedro_reducido, nodo_2)
                                    end


                                # SI TIENEN DOS CARAS EN COMUN
                                else

                                    # SI NINGUNA CARA COMUN ES DEL EXTERIOR DEL SOLIDO: se cuenta un cuarto
                                    if clasificacion_caras_hexaedro[caras_comunes] == [0,0]
                                        barra=conectividad_metamaterial[indice_1[i][1],3] 
                                        push!(barras_hexaedro_exacto, barra); push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                        push!(porcentajes_barras_dentro_hexaedro, 0.25)
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                    # SI LAS CARAS COMUN SON DEL EXTERIOR DEL SOLIDO: se cuenta entera
                                    elseif clasificacion_caras_hexaedro[caras_comunes] == [1,1]
                                        barra=conectividad_metamaterial[indice_1[i][1],3] 
                                        push!(barras_hexaedro_exacto, barra);push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                        push!(porcentajes_barras_dentro_hexaedro, 0.25) #1
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                    # SI DE LAS 2 CARAS COMUNES UNA CARA ES EXTERIOR Y LA OTRA NO: se cuenta la mitad
                                    else
                                        barra=conectividad_metamaterial[indice_1[i][1],3] 
                                        push!(barras_hexaedro_exacto, barra);push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                        push!(porcentajes_barras_dentro_hexaedro, 0.25) #0.5
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                    end
                                end


                            # SI EL NODO 2 ESTA EN EL EXTERIOR
                            else                             
                                coords_punto_corte, caras_pertenece_punto_corte=punto_corte_barras_caras_hexaedro(coords_nodo_1,coords_nodo_2,coords_vertices_hexaedro,volumen_hexaedro)

                                # SI LA BARRA QUE FORMAN CORTA A ALGUNA CARA: hay que determinar donde esta el punto de corte
                                if coords_punto_corte != 0

                                    # aumentamos las matrices de coordenadas (si el nodo no esta metido todavia) y de conectividad de metamaterial con el nodo del punto de corte
                                    # hay que ver si el nodo esta ya metido o no. Para ello, buscamos solo entre las coordenadas nuevas
                                    k=false; p=0
                                    for indice in N_nodos_metamaterial:length(coords_metamaterial[:,1]) 
                                        if abs(coords_punto_corte[1]-coords_metamaterial[indice,1])<epsilon
                                            if abs(coords_punto_corte[2]-coords_metamaterial[indice,2])<epsilon
                                                if abs(coords_punto_corte[3]-coords_metamaterial[indice,3])<epsilon
                                                    p=indice
                                                    k=true; break   
                                                else
                                                end
                                            else
                                            end
                                        else
                                        end
                                    end

                                    indice_conectiv_pto_corte=indice_conectiv_pto_corte+1
                                    # si no se ha registrado el nodo todavia, le registramos ahora
                                    if k==false
                                        indice_coords_pto_corte=length(coords_metamaterial[:,1])+1
                                        coords_metamaterial=vcat(coords_metamaterial,[coords_punto_corte indice_coords_pto_corte])
                                        conectividad_metamaterial=vcat(conectividad_metamaterial,[nodo_1 indice_coords_pto_corte indice_conectiv_pto_corte])
                                    # si ya se habia registrado el nodo
                                    else
                                        indice_coords_pto_corte=p
                                        conectividad_metamaterial=vcat(conectividad_metamaterial,[nodo_1 indice_coords_pto_corte indice_conectiv_pto_corte])
                                    end


                                    # Hay que determinar cuantas caras tienen en comun el punto de corte y el nodo 1 y cuales son
                                    caras_comunes=[] 
                                    for cara in caras_pertenece_nodo_1 
                                        if findall(in(cara), caras_pertenece_punto_corte) != []
                                            push!(caras_comunes, cara)
                                        else
                                        end
                                    end
                                    N_caras_comunes=length(caras_comunes)

                                    # SI NO TIENEN NINGUNA CARA EN COMUN EL PUNTO 1 Y EL PUNTO DE CORTE
                                    if N_caras_comunes==0
                                        push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                        push!(porcentajes_barras_dentro_hexaedro, 1)
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)  
                                        push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        #push!(Nodos_metam_hexaedro_reducido, nodo_2)


                                    # SI SOLO TIENEN UNA CARA EN COMUN
                                    elseif N_caras_comunes==1

                                        # SI LA CARA COMUN ES DEL EXTERIOR DEL SOLIDO: se cuenta entera
                                        if clasificacion_caras_hexaedro[caras_comunes] == [1]
                                            push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                            push!(porcentajes_barras_dentro_hexaedro, 0.5) #1
                                            push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)  
                                            push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                            #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                        # SI LA CARA COMUN NO ES DEL EXTERIOR DEL SOLIDO: se cuenta la mitad
                                        else
                                            push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                            push!(porcentajes_barras_dentro_hexaedro, 0.5)
                                            push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)  
                                            push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                            #push!(Nodos_metam_hexaedro_reducido, nodo_2)
                                        end


                                    # SI TIENEN DOS CARAS EN COMUN
                                    else

                                        # SI NINGUNA CARA COMUN ES DEL EXTERIOR DEL SOLIDO: se cuenta un cuarto
                                        if clasificacion_caras_hexaedro[caras_comunes] == [0,0]
                                            push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                            push!(porcentajes_barras_dentro_hexaedro, 0.25)
                                            push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)  
                                            push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                            #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                        # SI LAS CARAS COMUNES SON DEL EXTERIOR DEL SOLIDO: se cuenta entera
                                        elseif clasificacion_caras_hexaedro[caras_comunes] == [1,1]
                                            push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                            push!(porcentajes_barras_dentro_hexaedro, 0.25) #1
                                            push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)
                                            push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                            #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                        # SI DE LAS 2 CARAS COMUNES UNA CARA ES EXTERIOR Y LA OTRA NO: se cuenta la mitad
                                        else
                                            push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                            push!(porcentajes_barras_dentro_hexaedro, 0.25) #0.5
                                            push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1) 
                                            push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                                                 
                                            #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                        end
                                    end

                                # SI LA BARRA NO CORTA EN NINGUN PUNTO: no se cuenta
                                else
                                    push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1) 
                                     
                                end 

                            end

                        end
                    end
                end

            # SI NO ESTA EN EL HEXAEDRO: no se hace nada
            else
            end

        end
        
        unique!(sort!(Nodos_metam_hexaedro_exacto)); unique!(sort!(Nodos_metam_hexaedro_reducido))
        unique!(sort!(Nodos_metam_frontera_exacta_hexaedro)); unique!(sort!(Nodos_metam_frontera_extendido_hexaedro))
        Nodos_metam_hexaedro_exacto=unique!(sort!(vcat(Nodos_metam_hexaedro_exacto, Nodos_metam_frontera_exacta_hexaedro))); Nodos_metam_hexaedro_extendido=unique!(sort!(vcat(Nodos_metam_hexaedro_reducido, Nodos_metam_frontera_extendido_hexaedro)))
        indices_unicos_barras=unique(i->barras_hexaedro_exacto[i], 1:length(barras_hexaedro_exacto))

        for i in indices_unicos_barras
            barra=barras_hexaedro_exacto[i]
            porcentajes_barras_dentro_hexaedro_dict[barra]=porcentajes_barras_dentro_hexaedro[i]
        end
        barras_hexaedro_exacto=barras_hexaedro_exacto[indices_unicos_barras]; barras_hexaedro_extendido=barras_hexaedro_extendido[indices_unicos_barras] 
        unique!(sort!(barras_hexaedro_reducido)) 
        indices=findall(x-> x<=N_elementos_metam, barras_hexaedro_reducido)
        barras_hexaedro_reducido=barras_hexaedro_reducido[indices]


            # println(elemento_solido)
            # println("conectividad_ext")
            # println(barras_hexaedro_extendido);
            # println("porcentajes")
            # println(porcentajes_barras_dentro_hexaedro)
            # println("porcentajes dict")
            # println(porcentajes_barras_dentro_hexaedro_dict)
            # println("conectividad_int")
            # println(barras_hexaedro_reducido)
            # println("nodos_ext")
            # println(Nodos_metam_hexaedro_extendido)
            # println("nodos_int")
            # println(Nodos_metam_hexaedro_reducido)
            # println("frontera")
            # println(Nodos_metam_frontera_exacta_hexaedro)
            #println("conectividad")
            #println(conectividad_metamaterial[barras_hexaedro_extendido,:])
            #readline()

        # GENERAR LAS MATRICES DE CONECTIVIDAD (NORMAL Y EXTENDIDA) DE CADA ELEMENTO
        conectiv_metam_hexaedros_exacto[elemento_solido]=conectividad_metamaterial[barras_hexaedro_exacto,:]; conectiv_metam_hexaedros_extendido[elemento_solido]=conectividad_metamaterial[barras_hexaedro_extendido,:]; conectiv_metam_hexaedros_reducido[elemento_solido]=conectividad_metamaterial[barras_hexaedro_reducido,:]  # conectividad de elementos que tienen los dos nodos en el interior del hexaedro o en su frontera
        porcentajes_barras_dentro_hexaedros[elemento_solido]=porcentajes_barras_dentro_hexaedro_dict
        coords_metam_fronteras_exactas_hexaedros[elemento_solido]=coords_metamaterial[Nodos_metam_frontera_exacta_hexaedro,:]
        coords_metam_hexaedros_exacto[elemento_solido]=coords_metamaterial[Nodos_metam_hexaedro_exacto,:]; coords_metam_hexaedros_extendido[elemento_solido]=coords_metamaterial[Nodos_metam_hexaedro_extendido,:]; coords_metam_hexaedros_reducido[elemento_solido]=coords_metamaterial[Nodos_metam_hexaedro_reducido,:]

        # debemos calcular las tensiones de pandeo con los elementos completos (extendido), pero deben mantener la numeracion de las barras exactas y el mismo orden
        tensiones_pandeo_elementos[elemento_solido]=tensiones_criticas_pandeo(coords_metam_hexaedros_extendido[elemento_solido], [conectiv_metam_hexaedros_extendido[elemento_solido][:,1:2] conectiv_metam_hexaedros_exacto[elemento_solido][:,3]], Young_material(), radio)

        #println(conectiv_metam_hexaedros_extendido[elemento_solido]);readline()
    end

    return coords_metamaterial, conectividad_metamaterial, conectiv_metam_hexaedros_exacto, conectiv_metam_hexaedros_extendido, conectiv_metam_hexaedros_reducido, porcentajes_barras_dentro_hexaedros, coords_metam_hexaedros_exacto, coords_metam_hexaedros_extendido, coords_metam_hexaedros_reducido, coords_metam_fronteras_exactas_hexaedros, tensiones_pandeo_elementos

end #function nodos_metamaterial_en_elementos



# funcion para clasificar los nodos de metamaterial dentro de los hexaedros
# En este caso, las barras que se encuentran entre 2 elementos, se consideran solo en uno de ellos
function nodos_metamaterial_en_elementos_sin_repetir(X,Y,Z,coords_metamaterial,conectividad_metamaterial,coords_solido,conectividad_solido,volumenes_hexaedros,caras_hexaedros,radio)
    
    N_elementos_solidos=length(conectividad_solido[:,1]) # total de elementos hexaedricos
    N_nodos_metamaterial=length(coords_metamaterial[:,1]) # total de nodos en la geometria del metamaterial
    N_elementos_metam=length(conectividad_metamaterial[:,1])
    #indice_coords_pto_corte=N_nodos_metamaterial
    indice_conectiv_pto_corte=N_elementos_metam
    epsilon=1e-10; young_material=Young_material()
    barras_utilizadas=[]; elementos_barras_utlizadas=[]
    
    conectiv_metam_hexaedros_exacto=Dict(); conectiv_metam_hexaedros_extendido=Dict(); conectiv_metam_hexaedros_reducido=Dict()  # conectividad de elementos que tienen los dos nodos en el interior del hexaedro o en su frontera
    barras_ya_consideradas_en_otro_elemento_hexaedros=Dict();  # barras que ya se han considerado, el elemento en el que estan y elporcentaje de barra que le toca
    porcentajes_barras_dentro_hexaedros=Dict()      
    coords_metam_hexaedros_exacto=Dict(); coords_metam_hexaedros_extendido=Dict(); coords_metam_hexaedros_reducido=Dict()
    coords_metam_fronteras_exactas_hexaedros=Dict()
    caras_nodos_metam_exacto=Dict() # las caras en las que estan los nodos de la frontera
    tensiones_pandeo_elementos=Dict()
    registro_nodos_padres_hijos=Dict()


    for elemento_solido in 1:N_elementos_solidos

        barras_hexaedro_exacto=[];barras_hexaedro_extendido=[];barras_hexaedro_reducido=[]
        porcentajes_barras_dentro_hexaedro=[]; porcentajes_barras_dentro_hexaedro_dict=Dict()                      
        Nodos_metam_hexaedro_exacto=[];Nodos_metam_hexaedro_extendido=[];Nodos_metam_hexaedro_reducido=[] 
        Nodos_metam_frontera_exacta_hexaedro=[];Nodos_metam_frontera_extendido_hexaedro=[]    
        barras_ya_consideradas_en_otro_elemento=[]; hexaedro_barras_ya_consideradas_en_otro_elemento=[]
        porcentajes_dentro_hexaedro_barras_ya_consideradas_en_otro_elemento=[]
        hexaedro_barras_ya_consideradas_en_otro_elemento_unico=[]; porcentajes_dentro_hexaedro_barras_ya_consideradas_en_otro_elemento_unico=[]

        coords_vertices_hexaedro=zeros(8,3) #4 nodos con 3 coordenadas por nodo (vertices del hexaedro)
        for i in 1:8 
            nodo_solido=conectividad_solido[elemento_solido,i] # numero de nodo de los vertices del solido
            coords_vertices_hexaedro[i,:]=coords_solido[nodo_solido,1:3]  # coordenadas de los nodos del elemento hexaedrico
        end
        volumen_hexaedro = volumenes_hexaedros[elemento_solido] #volumen del hexaedro en el que estamos

        clasificacion_caras_hexaedro=tipo_caras_hexaedro(coords_vertices_hexaedro,X,Y,Z)
        #println(clasificacion_caras_hexaedro)


        # ESTUDIAR CADA NODO DE METAMATERIAL
        for nodo_1 in 1:N_nodos_metamaterial
            coords_nodo_1=coords_metamaterial[nodo_1,1:3]
            valor_1, caras_pertenece_nodo_1=is_in_hexaedro(coords_nodo_1, coords_vertices_hexaedro, volumen_hexaedro)


            # SI ESTA EN EL HEXAEDRO (INTERIOR O CARAS)
            if  valor_1 != 0 #para ver si esta dentro o no se usan y comparan los volumenes
                
                # buscar con que nodos tiene conectividad el nodo 1
                indice_1=findall(in(nodo_1),conectividad_metamaterial[:,1:2])
                N_nodos_contacto=length(indice_1) # numero total de conectividades del nodo

                # hay que asegurarse de que al menos alguna conectividad sea con otro nodo de dentro para que no este aislado
                nodos_2=zeros(Int64,N_nodos_contacto); valores_2=zeros(Int64,N_nodos_contacto); coords_nodos_2=zeros(N_nodos_contacto,3); caras_pertenecen_nodos_2=Dict()

                for i in 1:N_nodos_contacto
                    # obtenemos para cada nodo con el que tiene conectividad, cual es ese nodo
                    if indice_1[i][2]==1
                        nodos_2[i]=conectividad_metamaterial[indice_1[i][1],2]
                    else
                        nodos_2[i]=conectividad_metamaterial[indice_1[i][1],1]
                    end
                    coords_nodos_2[i,:]=coords_metamaterial[nodos_2[i],1:3]
                    valores_2[i], caras_pertenecen_nodos_2[i]=is_in_hexaedro(coords_nodos_2[i,:], coords_vertices_hexaedro, volumen_hexaedro)
                end
                

                push!(Nodos_metam_hexaedro_exacto, nodo_1); push!(Nodos_metam_hexaedro_reducido, nodo_1)

                # RECORREMOS TODOS LOS NODOS CON LOS QUE TIENE CONTACTO
                for i in 1:N_nodos_contacto

                    valor_2=valores_2[i]
                    nodo_2=nodos_2[i]
                    coords_nodo_2=coords_nodos_2[i,:]
                    caras_pertenece_nodo_2=caras_pertenecen_nodos_2[i]


                    # SI EL NODO 1 ESTA EN EL INTERIOR
                    if caras_pertenece_nodo_1 == [] # en este caso estará en el interior del hexaedro
                        
                        # SI EL NODO 2 TAMBIEN ESTA EN EL INTERIOR: se cuenta el elemento entero
                        if valor_2 == 1              
                            barra=conectividad_metamaterial[indice_1[i][1],3]    
                            push!(barras_hexaedro_exacto,barra); push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                            push!(porcentajes_barras_dentro_hexaedro, 1)

                        # SI EL NODO 2 ESTA EN LAS CARAS: se cuenta el elemento entero
                        elseif valor_2 == 2
                            barra=conectividad_metamaterial[indice_1[i][1],3] 
                            push!(barras_hexaedro_exacto, barra); push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                            push!(porcentajes_barras_dentro_hexaedro, 1)
                            push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)

                        # SI EL NODO 2 ESTA FUERA: se obtiene el punto de corte y se genera un nuevo elemento
                        else
                            # aumentamos las matrices de coordenadas (si el nuevo nodo no está ya metido) y de conectividad de metamaterial con el nodo del punto de corte
                            coords_punto_corte, caras_pertenece_punto_corte= punto_corte_barras_caras_hexaedro(coords_nodo_1,coords_nodo_2,coords_vertices_hexaedro,volumen_hexaedro)

                            # hay que ver si el nodo esta ya metido o no. Para ello, buscamos solo entre las coordenadas nuevas
                            k=false; p=0
                            for indice in N_nodos_metamaterial:length(coords_metamaterial[:,1]) 
                                if abs(coords_punto_corte[1]-coords_metamaterial[indice,1])<epsilon
                                    if abs(coords_punto_corte[2]-coords_metamaterial[indice,2])<epsilon
                                        if abs(coords_punto_corte[3]-coords_metamaterial[indice,3])<epsilon
                                            p=indice
                                            k=true; break   
                                        else
                                        end
                                    else
                                    end
                                else
                                end
                            end

                            indice_conectiv_pto_corte=length(conectividad_metamaterial[:,1])+1
                            # si no se ha registrado el nodo todavia, le registramos ahora
                            if k==false
                                indice_coords_pto_corte=length(coords_metamaterial[:,1])+1
                                coords_metamaterial=vcat(coords_metamaterial,[coords_punto_corte indice_coords_pto_corte])
                                conectividad_metamaterial=vcat(conectividad_metamaterial,[nodo_1 indice_coords_pto_corte indice_conectiv_pto_corte])
                            # si ya se habia registrado el nodo
                            else
                                indice_coords_pto_corte=p
                                conectividad_metamaterial=vcat(conectividad_metamaterial,[nodo_1 indice_coords_pto_corte indice_conectiv_pto_corte])
                            end
                            

                            push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                            push!(porcentajes_barras_dentro_hexaedro,1)
                            push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                            #push!(Nodos_metam_hexaedro_reducido, nodo_2)
                        

                            # Como el punto de corte esta en alguna cara de los hexaedros, obtenemos esa cara en numeracion global                           
                            caras_pto_corte_globales=caras_hexaedros[elemento_solido,caras_pertenece_punto_corte]
                            for cara in caras_pto_corte_globales
                                # hay que ver si ya hemos creado el vector para esa cara en el diccionario o no
                                if haskey(caras_nodos_metam_exacto, cara)==true
                                    push!(caras_nodos_metam_exacto[cara],indice_coords_pto_corte)
                                else
                                    caras_nodos_metam_exacto[cara]=[indice_coords_pto_corte]
                                end
                            end


                            # SI EL NODO 1 ESTA DENTRO Y EL 2 FUERA, PODRIA DARSE EL CASO DE QUE NO ESTUVIESEN EN ELEMENTOS ADYACENTES
                            # Tener un registro de a partir de que dos nodos originales (padres) se ha generado el intermedio (hijo)
                            # si no esta registrado, le registramos
                            if haskey(registro_nodos_padres_hijos,(nodo_2,nodo_1))==false
                                registro_nodos_padres_hijos[nodo_1,nodo_2]=indice_coords_pto_corte

                            # si ya esta registrado, hay que comprobar si lo esta con el mismo puto de corte u otro (de la misma barra)
                            else
                                # si es el mismo, no hacemos nada (la barra se divide solo entre 2 (lo normal))
                                if registro_nodos_padres_hijos[nodo_2,nodo_1]==indice_coords_pto_corte
                                    
                                # si no es el mismo, significa que la barra esta dividida entre 3 y hay que registrar el elemento intermedio extra que se genera entre los dos puntos de corte(que ademas estará en otro hexaedro)
                                else
                                    println("ERROR: barra en 3 elementos: Fallo de precision por solo considerar 2")
                                    indice_coords_pto_corte_2=registro_nodos_padres_hijos[nodo_2,nodo_1]
                                    indice_conectiv_pto_corte=length(conectividad_metamaterial[:,1])+1
                                    println("barra no considerada: ",indice_coords_pto_corte_2," ",indice_coords_pto_corte)
                                    #conectividad_metamaterial=vcat(conectividad_metamaterial,[indice_coords_pto_corte_2 indice_coords_pto_corte indice_conectiv_pto_corte])
                                end
                            end
                        end


                    # SI EL NODO 1 ESTA EN ALGUNA CARA
                    else 

                        # Si el nodo 1 esta en alguna cara de los hexaedros, obtenemos esa cara en numeracion global                           
                        caras_nodo_1_globales=caras_hexaedros[elemento_solido,caras_pertenece_nodo_1]
                        for cara in caras_nodo_1_globales
                            # hay que ver si ya hemos creado el vecor para esa cara en el diccionario o no
                            if haskey(caras_nodos_metam_exacto, cara)==true
                                push!(caras_nodos_metam_exacto[cara],nodo_1)
                            else
                                caras_nodos_metam_exacto[cara]=[nodo_1]
                            end
                        end


                        # SI EL NODO 2 ESTA EN EL INTERIOR: se cuenta el elemento entero
                        if valor_2 == 1
                            barra=conectividad_metamaterial[indice_1[i][1],3] 
                            push!(barras_hexaedro_exacto, barra); push!(barras_hexaedro_extendido, barra); push!(barras_hexaedro_reducido, barra)
                            push!(porcentajes_barras_dentro_hexaedro, 1)
                            push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)  
                            #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                        # SI EL NODO 2 ESTA EN LAS CARAS: se cuenta parte del elemento o todo entero
                        elseif valor_2 == 2

                            # Si el nodo 2 esta en alguna cara de los hexaedros, obtenemos esa cara en numeracion global                           
                            caras_nodo_2_globales=caras_hexaedros[elemento_solido,caras_pertenece_nodo_2]
                            for cara in caras_nodo_2_globales
                                # hay que ver si ya hemos creado el vecor para esa cara en el diccionario o no
                                if haskey(caras_nodos_metam_exacto, cara)==true
                                    push!(caras_nodos_metam_exacto[cara],nodo_2)
                                else
                                    caras_nodos_metam_exacto[cara]=[nodo_2]
                                end
                            end


                            # Hay que determinar cuantas caras tienen en comun y cuales son
                            caras_comunes=[] 
                            for cara in caras_pertenece_nodo_1 
                                if findall(in(cara), caras_pertenece_nodo_2) != []
                                    push!(caras_comunes, cara)
                                else
                                end
                            end
                            N_caras_comunes=length(caras_comunes)

                            # SI NO TIENEN NINGUNA CARA EN COMUN
                            if N_caras_comunes==0
                                barra=conectividad_metamaterial[indice_1[i][1],3] 
                                push!(barras_hexaedro_exacto, barra);push!(barras_hexaedro_extendido, barra); push!(barras_hexaedro_reducido, barra)
                                push!(porcentajes_barras_dentro_hexaedro, 1)
                                push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                #push!(Nodos_metam_hexaedro_reducido, nodo_2)


                            # SI SOLO TIENEN UNA CARA EN COMUN
                            elseif N_caras_comunes==1

                                # SI LA CARA COMUN ES DEL EXTERIOR DEL SOLIDO: se cuenta entero el elemento
                                if clasificacion_caras_hexaedro[caras_comunes] == [1]
                                    barra=conectividad_metamaterial[indice_1[i][1],3] 
                                    push!(barras_hexaedro_exacto, barra); push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                    push!(porcentajes_barras_dentro_hexaedro, 0.5) #1
                                    push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                    push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                    #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                # SI LA CARA COMUN NO ES DEL EXTERIOR DEL SOLIDO: se cuenta la mitad
                                else
                                    barra=conectividad_metamaterial[indice_1[i][1],3]
                                    if (indice=findall(in(barra),barras_utilizadas))==[]  
                                        push!(barras_hexaedro_exacto, barra); push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                        push!(barras_utilizadas,barra); push!(elementos_barras_utlizadas,elemento_solido)
                                        push!(porcentajes_barras_dentro_hexaedro, 0.5)
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        
                                    elseif findall(in(barra),barras_hexaedro_exacto)==[]  
                                        push!(barras_ya_consideradas_en_otro_elemento, barra)
                                        push!(hexaedro_barras_ya_consideradas_en_otro_elemento, elementos_barras_utlizadas[indice][1])
                                        push!(porcentajes_dentro_hexaedro_barras_ya_consideradas_en_otro_elemento, 0.5)
                                    else
                                    end                     
                                    #push!(Nodos_metam_hexaedro_reducido, nodo_2)
                                end


                            # SI TIENEN DOS CARAS EN COMUN
                            else

                                # SI NINGUNA CARA COMUN ES DEL EXTERIOR DEL SOLIDO: se cuenta un cuarto
                                if clasificacion_caras_hexaedro[caras_comunes] == [0,0]
                                    
                                    barra=conectividad_metamaterial[indice_1[i][1],3] 
                                    if (indice=findall(in(barra),barras_utilizadas))==[]  
                                        push!(barras_hexaedro_exacto, barra); push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                        push!(barras_utilizadas,barra); push!(elementos_barras_utlizadas,elemento_solido)
                                        push!(porcentajes_barras_dentro_hexaedro, 0.25)
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        #push!(Nodos_metam_hexaedro_reducido, nodo_2)
                                    elseif findall(in(barra),barras_hexaedro_exacto)==[] 
                                        push!(barras_ya_consideradas_en_otro_elemento, barra)
                                        push!(hexaedro_barras_ya_consideradas_en_otro_elemento, elementos_barras_utlizadas[indice][1])
                                        push!(porcentajes_dentro_hexaedro_barras_ya_consideradas_en_otro_elemento, 0.25)
                                    else
                                    end   
                                    

                                # SI LAS CARAS COMUN SON DEL EXTERIOR DEL SOLIDO: se cuenta entera
                                elseif clasificacion_caras_hexaedro[caras_comunes] == [1,1]
                                    barra=conectividad_metamaterial[indice_1[i][1],3] 
                                    push!(barras_hexaedro_exacto, barra);push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                    push!(porcentajes_barras_dentro_hexaedro, 0.25) #1
                                    push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                    push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                    #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                # SI DE LAS 2 CARAS COMUNES UNA CARA ES EXTERIOR Y LA OTRA NO: se cuenta la mitad
                                else
                                    barra=conectividad_metamaterial[indice_1[i][1],3]   
                                    
                                    if (indice=findall(in(barra),barras_utilizadas))==[]  
                                        push!(barras_hexaedro_exacto, barra); push!(barras_hexaedro_extendido, barra) ;push!(barras_hexaedro_reducido, barra)
                                        push!(barras_utilizadas,barra); push!(elementos_barras_utlizadas,elemento_solido)
                                        push!(porcentajes_barras_dentro_hexaedro, 0.25) #0.5
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_2); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        #push!(Nodos_metam_hexaedro_reducido, nodo_2)
                                    elseif findall(in(barra),barras_hexaedro_exacto)==[] 
                                        push!(barras_ya_consideradas_en_otro_elemento, barra)
                                        push!(hexaedro_barras_ya_consideradas_en_otro_elemento, elementos_barras_utlizadas[indice][1])
                                        push!(porcentajes_dentro_hexaedro_barras_ya_consideradas_en_otro_elemento, 0.25) #0.5
                                    else
                                    end   
                                end
                            end


                        # SI EL NODO 2 ESTA EN EL EXTERIOR
                        else                             
                            coords_punto_corte, caras_pertenece_punto_corte=punto_corte_barras_caras_hexaedro(coords_nodo_1,coords_nodo_2,coords_vertices_hexaedro,volumen_hexaedro)

                            # SI LA BARRA QUE FORMAN CORTA A ALGUNA CARA: hay que determinar donde esta el punto de corte
                            if coords_punto_corte != 0

                                # aumentamos las matrices de coordenadas (si el nodo no esta metido todavia) y de conectividad de metamaterial con el nodo del punto de corte
                                # hay que ver si el nodo esta ya metido o no. Para ello, buscamos solo entre las coordenadas nuevas
                                k=false; p=0
                                for indice in N_nodos_metamaterial:length(coords_metamaterial[:,1]) 
                                    if abs(coords_punto_corte[1]-coords_metamaterial[indice,1])<epsilon
                                        if abs(coords_punto_corte[2]-coords_metamaterial[indice,2])<epsilon
                                            if abs(coords_punto_corte[3]-coords_metamaterial[indice,3])<epsilon
                                                p=indice
                                                k=true; break   
                                            else
                                            end
                                        else
                                        end
                                    else
                                    end
                                end


                                salida=false; element=0
                                # si no se ha registrado el nodo todavia, le registramos ahora
                                if k==false
                                    indice_conectiv_pto_corte=length(conectividad_metamaterial[:,1])+1
                                    indice_coords_pto_corte=length(coords_metamaterial[:,1])+1
                                    coords_metamaterial=vcat(coords_metamaterial,[coords_punto_corte indice_coords_pto_corte])
                                    conectividad_metamaterial=vcat(conectividad_metamaterial,[nodo_1 indice_coords_pto_corte indice_conectiv_pto_corte])
                                
                                # si ya se habia registrado el nodo
                                else
                                    indice_coords_pto_corte=p
                                    
                                    # vemos si el elemento ya esta registrado o no                                        
                                    for indice_conectiv in N_elementos_metam:length(conectividad_metamaterial[:,1]) 
                                        if nodo_1==conectividad_metamaterial[indice_conectiv,1]
                                            if indice_coords_pto_corte==conectividad_metamaterial[indice_conectiv,2]
                                                # si ya estaba registrado, no hacemos nada
                                                salida==true; element=indice_conectiv
                                                break
                                            else
                                            end
                                        else
                                        end
                                    end

                                    # si ya estaba registrado, no hacemos nada
                                    if salida==true
                                        indice_conectiv_pto_corte=element
                                    # si no estaba registrado, lo registramos
                                    else
                                        indice_conectiv_pto_corte=length(conectividad_metamaterial[:,1])+1
                                        conectividad_metamaterial=vcat(conectividad_metamaterial,[nodo_1 indice_coords_pto_corte indice_conectiv_pto_corte])
                                    end
                                end


                                
                                # Como el punto de corte esta en alguna cara de los hexaedros, obtenemos esa cara en numeracion global                           
                                caras_pto_corte_globales=caras_hexaedros[elemento_solido,caras_pertenece_punto_corte]
                                for cara in caras_pto_corte_globales
                                    # hay que ver si ya hemos creado el vector para esa cara en el diccionario o no
                                    if haskey(caras_nodos_metam_exacto, cara)==true
                                        push!(caras_nodos_metam_exacto[cara],indice_coords_pto_corte)
                                    else
                                        caras_nodos_metam_exacto[cara]=[indice_coords_pto_corte]
                                    end
                                end


                                # Hay que determinar cuantas caras tienen en comun el punto de corte y el nodo 1 y cuales son
                                caras_comunes=[] 
                                for cara in caras_pertenece_nodo_1 
                                    if findall(in(cara), caras_pertenece_punto_corte) != []
                                        push!(caras_comunes, cara)
                                    else
                                    end
                                end
                                N_caras_comunes=length(caras_comunes)


                                # SI NO TIENEN NINGUNA CARA EN COMUN EL PUNTO 1 Y EL PUNTO DE CORTE
                                if N_caras_comunes==0
                                    push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                    push!(porcentajes_barras_dentro_hexaedro, 1)
                                    push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)  
                                    push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                    #push!(Nodos_metam_hexaedro_reducido, nodo_2)


                                # SI SOLO TIENEN UNA CARA EN COMUN
                                elseif N_caras_comunes==1

                                    # SI LA CARA COMUN ES DEL EXTERIOR DEL SOLIDO: se cuenta entera
                                    if clasificacion_caras_hexaedro[caras_comunes] == [1]
                                        push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                        push!(porcentajes_barras_dentro_hexaedro, 0.5) #1
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)  
                                        push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                    # SI LA CARA COMUN NO ES DEL EXTERIOR DEL SOLIDO: se cuenta la mitad
                                    else
                                        # si el elemento se ha creado nuevo
                                        barra=indice_conectiv_pto_corte
                                        if (indice=findall(in(barra),barras_utilizadas))==[]   
                                            push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                            push!(barras_utilizadas,indice_conectiv_pto_corte); push!(elementos_barras_utlizadas,elemento_solido)
                                            push!(porcentajes_barras_dentro_hexaedro, 0.5)
                                            push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)  
                                            push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2) 
                                            #push!(Nodos_metam_hexaedro_reducido, nodo_2)
                                        else
                                            push!(barras_ya_consideradas_en_otro_elemento, indice_conectiv_pto_corte)
                                            push!(hexaedro_barras_ya_consideradas_en_otro_elemento, elementos_barras_utlizadas[indice][1])
                                            push!(porcentajes_dentro_hexaedro_barras_ya_consideradas_en_otro_elemento, 0.5)
                                        end
                                        
                                    end


                                # SI TIENEN DOS CARAS EN COMUN
                                else

                                    # SI NINGUNA CARA COMUN ES DEL EXTERIOR DEL SOLIDO: se cuenta un cuarto
                                    if clasificacion_caras_hexaedro[caras_comunes] == [0,0]
                                        # si el elemento se ha creado nuevo
                                        barra=indice_conectiv_pto_corte
                                        if (indice=findall(in(barra),barras_utilizadas))==[] 
                                            push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                            push!(barras_utilizadas,indice_conectiv_pto_corte); push!(elementos_barras_utlizadas,elemento_solido)
                                            push!(porcentajes_barras_dentro_hexaedro, 0.25)
                                            push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)  
                                            push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2) 
                                        else
                                            push!(barras_ya_consideradas_en_otro_elemento, indice_conectiv_pto_corte)
                                            push!(hexaedro_barras_ya_consideradas_en_otro_elemento, elementos_barras_utlizadas[indice][1])
                                            push!(porcentajes_dentro_hexaedro_barras_ya_consideradas_en_otro_elemento, 0.25) #0.5
                                        end

                                    # SI LAS CARAS COMUN SON DEL EXTERIOR DEL SOLIDO: se cuenta entera
                                    elseif clasificacion_caras_hexaedro[caras_comunes] == [1,1]
                                        push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                        push!(porcentajes_barras_dentro_hexaedro, 0.25) #1
                                        push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)
                                        push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2)  
                                        #push!(Nodos_metam_hexaedro_reducido, nodo_2)

                                    # SI DE LAS 2 CARAS COMUNES UNA CARA ES EXTERIOR Y LA OTRA NO: se cuenta la mitad
                                    else
                                        # si el elemento se ha creado nuevo
                                        barra=indice_conectiv_pto_corte
                                        if (indice=findall(in(barra),barras_utilizadas))==[] 
                                            push!(barras_hexaedro_exacto, indice_conectiv_pto_corte); push!(barras_hexaedro_extendido, conectividad_metamaterial[indice_1[i][1],3])
                                            push!(barras_utilizadas,indice_conectiv_pto_corte); push!(elementos_barras_utlizadas,elemento_solido)
                                            push!(porcentajes_barras_dentro_hexaedro, 0.25) #0.5
                                            push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1)  
                                            push!(Nodos_metam_frontera_exacta_hexaedro, indice_coords_pto_corte); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_2) 
                                            #push!(Nodos_metam_hexaedro_reducido, nodo_2)
                                        else
                                            push!(barras_ya_consideradas_en_otro_elemento, indice_conectiv_pto_corte)
                                            push!(hexaedro_barras_ya_consideradas_en_otro_elemento, elementos_barras_utlizadas[indice][1])
                                            push!(porcentajes_dentro_hexaedro_barras_ya_consideradas_en_otro_elemento, 0.25) #0.5
                                        end
                                        

                                    end
                                end

                            # SI LA BARRA NO CORTA EN NINGUN PUNTO: no se cuenta
                            else
                                push!(Nodos_metam_frontera_exacta_hexaedro, nodo_1); push!(Nodos_metam_frontera_extendido_hexaedro, nodo_1) 
                                    
                            end 

                        end

                    end
                end

            # SI NO ESTA EN EL HEXAEDRO: no se hace nada
            else
            end

        end
        
        unique!(sort!(Nodos_metam_hexaedro_exacto)); unique!(sort!(Nodos_metam_hexaedro_reducido))
        unique!(sort!(Nodos_metam_frontera_exacta_hexaedro)); unique!(sort!(Nodos_metam_frontera_extendido_hexaedro))
        Nodos_metam_hexaedro_exacto=unique!(sort!(vcat(Nodos_metam_hexaedro_exacto, Nodos_metam_frontera_exacta_hexaedro))); Nodos_metam_hexaedro_extendido=unique!(sort!(vcat(Nodos_metam_hexaedro_reducido, Nodos_metam_frontera_extendido_hexaedro)))
        indices_unicos_barras=unique(i->barras_hexaedro_exacto[i], 1:length(barras_hexaedro_exacto))

        # Nos quedamos con las barras sin repetir
        for i in indices_unicos_barras
            barra=barras_hexaedro_exacto[i]
            porcentajes_barras_dentro_hexaedro_dict[barra]=porcentajes_barras_dentro_hexaedro[i]
        end
        barras_hexaedro_exacto=barras_hexaedro_exacto[indices_unicos_barras]; barras_hexaedro_extendido=barras_hexaedro_extendido[indices_unicos_barras] 
        unique!(sort!(barras_hexaedro_reducido)) 
        indices=findall(x-> x<=N_elementos_metam, barras_hexaedro_reducido)
        barras_hexaedro_reducido=barras_hexaedro_reducido[indices] 

        indices_unicos_barras_ya_consideradas=unique(i->barras_ya_consideradas_en_otro_elemento[i], 1:length(barras_ya_consideradas_en_otro_elemento))
        coords_nodo_1=zeros(length(indices_unicos_barras_ya_consideradas),3); coords_nodo_2=zeros(length(indices_unicos_barras_ya_consideradas),3) 
        k=1
        for i in indices_unicos_barras_ya_consideradas
            barra= barras_ya_consideradas_en_otro_elemento[i]
            push!(porcentajes_dentro_hexaedro_barras_ya_consideradas_en_otro_elemento_unico, porcentajes_dentro_hexaedro_barras_ya_consideradas_en_otro_elemento[i])
            push!(hexaedro_barras_ya_consideradas_en_otro_elemento_unico, hexaedro_barras_ya_consideradas_en_otro_elemento[i])
            
            nodo_1=conectividad_metamaterial[barra,1]; coords_nodo_1[k,:]=coords_metamaterial[nodo_1, 1:3]'
            nodo_2=conectividad_metamaterial[barra,2]; coords_nodo_2[k,:]=coords_metamaterial[nodo_2, 1:3]'
            k=k+1
        end
        barras_ya_consideradas_en_otro_elemento=barras_ya_consideradas_en_otro_elemento[indices_unicos_barras_ya_consideradas]
        # generamos una unica matriz donde estan las barras ya utilizadas, sus porcentajes y en que hexaedro estan
        barras_ya_consideradas_en_otro_elemento=hcat(barras_ya_consideradas_en_otro_elemento, hexaedro_barras_ya_consideradas_en_otro_elemento_unico, porcentajes_dentro_hexaedro_barras_ya_consideradas_en_otro_elemento_unico, coords_nodo_1, coords_nodo_2)
        
        #println(elemento_solido)
        #println(barras_hexaedro_exacto);println()
        #println(barras_ya_consideradas_en_otro_elemento);readline()

            # println(elemento_solido)
            # println("conectividad_ext")
            # println(barras_hexaedro_extendido);
            # println("porcentajes")
            # println(porcentajes_barras_dentro_hexaedro)
            # println("porcentajes dict")
            # println(porcentajes_barras_dentro_hexaedro_dict)
            # println("conectividad_int")
            # println(barras_hexaedro_reducido)
            # println("nodos_ext")
            # println(Nodos_metam_hexaedro_extendido)
            # println("nodos_int")
            # println(Nodos_metam_hexaedro_reducido)
            # println("frontera")
            # println(Nodos_metam_frontera_exacta_hexaedro)
            #println("conectividad")
            #println(conectividad_metamaterial[barras_hexaedro_extendido,:])
            #readline()

        # GENERAR LAS MATRICES DE CONECTIVIDAD (NORMAL Y EXTENDIDA) DE CADA ELEMENTO
        conectiv_metam_hexaedros_exacto[elemento_solido]=conectividad_metamaterial[barras_hexaedro_exacto,:]; conectiv_metam_hexaedros_extendido[elemento_solido]=conectividad_metamaterial[barras_hexaedro_extendido,:]; conectiv_metam_hexaedros_reducido[elemento_solido]=conectividad_metamaterial[barras_hexaedro_reducido,:]  # conectividad de elementos que tienen los dos nodos en el interior del hexaedro o en su frontera
        porcentajes_barras_dentro_hexaedros[elemento_solido]=porcentajes_barras_dentro_hexaedro_dict
        barras_ya_consideradas_en_otro_elemento_hexaedros[elemento_solido]=barras_ya_consideradas_en_otro_elemento
        coords_metam_fronteras_exactas_hexaedros[elemento_solido]=coords_metamaterial[Nodos_metam_frontera_exacta_hexaedro,:]
        coords_metam_hexaedros_exacto[elemento_solido]=coords_metamaterial[Nodos_metam_hexaedro_exacto,:]; coords_metam_hexaedros_extendido[elemento_solido]=coords_metamaterial[Nodos_metam_hexaedro_extendido,:]; coords_metam_hexaedros_reducido[elemento_solido]=coords_metamaterial[Nodos_metam_hexaedro_reducido,:]

        # debemos calcular las tensiones de pandeo con los elementos completos (extendido), pero deben mantener la numeracion de las barras exactas y el mismo orden
        tensiones_pandeo_elementos[elemento_solido]=tensiones_criticas_pandeo(coords_metam_hexaedros_extendido[elemento_solido], [conectiv_metam_hexaedros_extendido[elemento_solido][:,1:2] conectiv_metam_hexaedros_exacto[elemento_solido][:,3]], Young_material(), radio)

        #println(conectiv_metam_hexaedros_extendido[elemento_solido])
        #readline()
    end

    # Quitamos los valores repetidos de las caras a las que pertenecen los nodos
    for cara in 1:length(caras_nodos_metam_exacto)
        caras_nodos_metam_exacto[cara]=unique!(sort!(caras_nodos_metam_exacto[cara]))
        #println(caras_nodos_metam_exacto[cara])
    end

    return coords_metamaterial, conectividad_metamaterial, conectiv_metam_hexaedros_exacto, conectiv_metam_hexaedros_extendido, conectiv_metam_hexaedros_reducido, porcentajes_barras_dentro_hexaedros, coords_metam_hexaedros_exacto, coords_metam_hexaedros_extendido, coords_metam_hexaedros_reducido, coords_metam_fronteras_exactas_hexaedros, caras_nodos_metam_exacto, barras_ya_consideradas_en_otro_elemento_hexaedros, tensiones_pandeo_elementos

end #function nodos_metamaterial_en_elementos


# ver que nodos de metamaterial estan en cada hexaedro
function clasificar_nodos_metam_en_hexaedros(coords_metam, coords_3D, conectividad_3D, volumenes_hexaedros)

    nodos_metam_en_hexaedros=Dict()

    for elemento_solido in 1:length(conectividad_3D[:,3])
        nodos_metam_en_hexaedros[elemento_solido]=[]

        coords_vertices_hexaedro=zeros(8,3)
        for i in 1:8 
            nodo_solido=conectividad_3D[elemento_solido,i] # numero de nodo de los vertices del solido
            coords_vertices_hexaedro[i,:]=coords_3D[nodo_solido,1:3]  # coordenadas de los nodos del elemento hexaedrico
        end

        volumen_hexaedro=volumenes_hexaedros[elemento_solido]

        for nodo in 1:length(coords_metam[:,4])

            coords_nodo=coords_metam[nodo,1:3]
            valor, caras_pertenece_nodo = is_in_hexaedro(coords_nodo, coords_vertices_hexaedro, volumen_hexaedro)

            if valor!=0
                push!(nodos_metam_en_hexaedros[elemento_solido], nodo)
            else
            end
        end

    end

    return nodos_metam_en_hexaedros

end



function punto_corte_barras_caras_hexaedro(coords_nodo1,coords_nodo2,coords_vertices_hexaedro,volumen_hexaedro) #esta funcion permite determinar si un elemento viga corta a alguna de las 4 caras del hexaedro
    # esta funcion es necesaria para los elementos viga que salen de un nodo que esta en una cara del hexaedro y que terminan en un nodo exterior...
    # ... porque podria ser que ese elemento atravesase al hexaedro y por eso acabase en un nodo de fuera
    
    # cualquier punto de una recta: [x,y,z]=[xa,ya,za]+lambda*[xb-xa,yb-ya,zb-za]
    # cualquier punto de un plano:  [x,y,z]=[x1,y1,z1]+alpha*[x2-x1,y2-y1,z2-z1]+beta*[x3-x1,y3-y1,z3-z1]
    # hay que determinar alpha,beta,lambda para determinar el punto de corte
    punto_corte=0; caras_pertenece_punto_corte=0

    xa,ya,za=coords_nodo1
    xb,yb,zb=coords_nodo2

    aux=[1 2 4 3; 1 5 6 2; 3 4 8 7; 1 3 7 5; 2 6 8 4; 5 7 8 6] # este criterio de eleccion de los vertices debe coincidir con el de las areas de las caras, para que se puedan comparar los valores
    for m in 1:6    
        x1,y1,z1=coords_vertices_hexaedro[aux[m,1],:]
        x2,y2,z2=coords_vertices_hexaedro[aux[m,2],:]
        x3,y3,z3=coords_vertices_hexaedro[aux[m,3],:]
        coords_vertices_cara_hexaedro=[x1 y1 z1; x2 y2 z2; x3 y3 z3]

        A=SMatrix{3, 3, Float64,9}([x2-x1  x3-x1  -(xb-xa);
                                    y2-y1  y3-y1  -(yb-ya);
                                    z2-z1  z3-z1  -(zb-za)])  #matriz creada con StaticArrays para ganar velocidad en la inversa
        C=[xa-x1, ya-y1, za-z1]
        alpha,beta,lambda= inv(A)*C

        epsilon=1e-12 # Para evitar el error numerico
        if (lambda > epsilon) && (lambda < (1-epsilon))  # en este caso, el punto de corte estaria fuera de la barra o justo en su extremo, asi que no lo consideramos
            punto_corte= Vector{Any}([xa,ya,za]+lambda*[xb-xa,yb-ya,zb-za]) # en este caso el punto de corte esta dentro del elemento barra
            #println("corte= ",punto_corte)

            ### ahora hay que ver si ese punto esta en la cara del hexaedro o se ha salido de la cara (aunque este en el mismo plano)
            valor_punto_corte, caras_pertenece_punto_corte = is_in_hexaedro(punto_corte,coords_vertices_hexaedro,volumen_hexaedro)
            if valor_punto_corte == 2
                punto_corte= punto_corte
                #punto_corte=lambda  # con que se cumpla para alguna de las caras, se considerará el elemento
                break
            else
                punto_corte=0
            end

        else
        end
    end

    #println(interseca_en_cara)
    #println("ver si interseca ",coords_nodo1," ",coords_nodo2," ",valor);readline()

    return transpose(punto_corte), caras_pertenece_punto_corte

end # function barra_interseca_caras_hexaedro



function clasificacion_gdl_multiescala(X,Y,Z,N_3D_x,N_3D_y,N_3D_z, coords_metam, conectiv_metam_hexaedros_exacto, coords_metam_hexaedros_exacto, coords_metam_fronteras_exactas_hexaedros, punto_restriccion, gdl_restriccion, punto_desplazamiento, inc_valor_desplazamiento)    

    N_elementos_solidos=length(coords_metam_hexaedros_exacto)

    # Pasar las restricciones del solido a restricciones del metamaterial (con 6 gdl en vez de 3)
    N_restricciones = length(gdl_restriccion[:,1])
    N_desplazamientos = length(inc_valor_desplazamiento[:,1])

    punto_restriccion_metam = punto_restriccion
    gdl_restriccion_metam = [gdl_restriccion zeros(N_restricciones,3)]
    punto_desplazamiento_metam = punto_desplazamiento
    inc_valor_desplazamiento_metam = [inc_valor_desplazamiento zeros(N_desplazamientos,3)]


    # GRADOS DE LIBERTAD DEL CONTORNO DEL SOLIDO COMPLETO QUE NO ESTEN EN LAS FRONTERAS ENTRE HEXAEDROS
    Nodos_metam_contorno_solido, gdl_contorno_solido=gdl_metamaterial_contorno_solido_no_fronteras(X,Y,Z,N_3D_x,N_3D_y,N_3D_z,coords_metam)
    # for i in Nodos_metam_contorno_solido
    #     println(coords_metam[i,:])
    # end

    # GRADOS DE LIBERTAD CON ALGUNA RESTRICCION Y EL VALOR DE ESA RESTRICCION (PUEDEN SER DEL CONTORNO O NO): SERAN LOS GDL FIJOS
    gdl_restring_globales, inc_valor_restricciones_fijos_global=condic_contorno(punto_restriccion_metam, gdl_restriccion_metam ,punto_desplazamiento_metam ,inc_valor_desplazamiento_metam ,coords_metam)
    #println(gdl_restring_globales, inc_valor_restricciones_fijos_global)


    gdl_zona_hexaedros=Dict()
    gdl_zona_todas_restricciones_hexaedros=Dict()
    gdl_zona_restricciones_fijas_hexaedros=Dict()
    inc_valor_restricciones_fijas_hexaedro=Dict()
    gdl_zona_restricciones_no_fijas_hexaedros=Dict()
    gdl_zona_libres_hexaedros=Dict()
    gdl_restricciones_no_fijas_global=[]


    for elemento_solido in 1:N_elementos_solidos
        
        # TODOS LOS GDL DE LA ZONA
        gdl_zona=[]
        for nodo in coords_metam_hexaedros_exacto[elemento_solido][:,4]
            append!(gdl_zona, collect(6*nodo-5:6*nodo))
        end

        # GRADOS DE LIBERTAD CON RESTRICCIONES FIJAS (NUNCA SE CAMBIAN) (PUEDEN SER DE LA FRONTERA O NO)
        indices=findall(in(gdl_zona),gdl_restring_globales)
        gdl_zona_restricciones_fijas=gdl_restring_globales[indices]
        inc_valor_restricciones_fijas_zona=inc_valor_restricciones_fijos_global[indices]


        # GRADOS DE LIBERTAD CON RESTRICCIONES NO FIJAS (SON LOS DE LA FRONTERA PERO SIN RESTRICCIONES FIJAS NI ESTAR EN EL CONTORNO EXTERIOR DEL SOLIDO) (SON LOS QUE SE VAN AJUSTANDO HASTA LA CONVERGENCIA DE DESPLAZAMIENTOS)
        gdl_frontera=[]
        for nodo in coords_metam_fronteras_exactas_hexaedros[elemento_solido][:,4]
            append!(gdl_frontera, collect(6*nodo-5:6*nodo))
        end

        gdl_zona_restricciones_no_fijas=[]
        for gdl in gdl_frontera
            if findall(in(gdl), gdl_contorno_solido)==[]
                if findall(in(gdl), gdl_zona_restricciones_fijas)==[]
                    push!(gdl_zona_restricciones_no_fijas, gdl)
                else
                end
            else
            end
        end


       
        # GRADOS DE LIBERTAD LIBRES
        gdl_zona_libres=[]
        for gdl in gdl_zona
            if findall(in(gdl), gdl_zona_restricciones_fijas)==[]
                if findall(in(gdl), gdl_zona_restricciones_no_fijas)==[]
                    push!(gdl_zona_libres, gdl)
                else
                end
            else
            end
        end        


        # GRADOS DE LIBERTAD CORRESPONDIENTES A TODAS LAS RESTRICCIONES (FIJAS Y NO FIJAS) (TODOS LOS QUE NO SON LIBRES)
        indices=findall(!in(gdl_zona_libres), gdl_zona)
        gdl_zona_restricciones_total=gdl_zona[indices]


        ### PASARLO A DICCIONARIOS CON TODOS LOS ELEMENTOS SOLIDOS 
        gdl_zona_hexaedros[elemento_solido]=gdl_zona
        gdl_zona_todas_restricciones_hexaedros[elemento_solido]=gdl_zona_restricciones_total
        gdl_zona_libres_hexaedros[elemento_solido]=gdl_zona_libres
        gdl_zona_restricciones_fijas_hexaedros[elemento_solido]=gdl_zona_restricciones_fijas
        inc_valor_restricciones_fijas_hexaedro[elemento_solido]=inc_valor_restricciones_fijas_zona
        gdl_zona_restricciones_no_fijas_hexaedros[elemento_solido]=gdl_zona_restricciones_no_fijas
        
        append!(gdl_restricciones_no_fijas_global, gdl_zona_restricciones_no_fijas)

    end

    gdl_restricciones_no_fijas_global=unique!(sort!(gdl_restricciones_no_fijas_global))

    return gdl_zona_hexaedros, gdl_zona_todas_restricciones_hexaedros, gdl_restricciones_no_fijas_global, gdl_zona_restricciones_fijas_hexaedros, inc_valor_restricciones_fijas_hexaedro, gdl_zona_restricciones_no_fijas_hexaedros, gdl_zona_libres_hexaedros

end



function clasificacion_gdl_zona_relajacion(elemento_solido, X,Y,Z,N_3D_x,N_3D_y,N_3D_z, coords_metam, N_capas_alrededor, hexaedros_alrededor_1_capa, caras_hexaedros, caras_contorno_solido, caras_nodos_metam_exacto, coords_metam_hexaedros_exacto, coords_metam_fronteras_exactas_hexaedros, gdl_restring_globales, gdl_contorno_solido_completo_no_fronteras)
    

    ### OBTENER LOS ELEMENTOS QUE ESTAN ALREDEDOR DEL QUE ESTUDIAMOS (TENIENDO EN CUENTA EL NUMERO DE CAPAS): INCLUIMOS AL PROPIO ELEMENTO
    elementos_alrededor=[elemento_solido]
    if N_capas_alrededor=="todas"
        N_capas_alrededor=maximum([N_3D_x, N_3D_y, N_3D_z])
    else
    end

    for capa in 1:N_capas_alrededor
        elementos_alrededor_nuevo=[]

        for elemento in elementos_alrededor
            #println(elemento)
            #println(hexaedros_alrededor_1_capa[elemento])
            append!(elementos_alrededor_nuevo, hexaedros_alrededor_1_capa[elemento])
            #readline()
        end

        append!(elementos_alrededor, elementos_alrededor_nuevo)
        unique!(sort!(elementos_alrededor))
    end


    ### OBTENER LOS GDL ASOCIADOS AL CONTORNO DE ESTA ZONA
    caras_totales_zona=[]; caras_repetidas_zona=[]; caras_contorno_zona=[]

    # OBTENER LAS CARAS QUE SON DEL CONTORNO DE LA ZONA (SERAN LAS QUE NO SE REPITAN)
    for elemento in elementos_alrededor
        
        for cara in caras_hexaedros[elemento,:]
            if findall(in(cara), caras_totales_zona)==[]
                append!(caras_totales_zona, cara)
            else
                append!(caras_repetidas_zona, cara)
            end
        end

        indices_caras_no_repetidas=findall(!in(caras_repetidas_zona), caras_totales_zona)
        caras_contorno_zona=caras_totales_zona[indices_caras_no_repetidas]
    end

    # OBTENER LOS NODOS Y GDL DE ESA FRONTERA
    gdl_contorno_zona=[]

    for cara in caras_contorno_zona
        for nodo in caras_nodos_metam_exacto[cara]
            append!(gdl_contorno_zona, collect(6*nodo-5:6*nodo))
        end
    end
    unique!(sort!(gdl_contorno_zona))

    # OBTENER LOS GDL CORRESPONDIENTES A LAS CARAS QUE SE REPITEN
    gdl_caras_repetidas=[]
    for cara in caras_repetidas_zona
        for nodo in caras_nodos_metam_exacto[cara]
            append!(gdl_caras_repetidas, collect(6*nodo-5:6*nodo))
        end
    end


    ### OBTENER TODOS LOS GDL DE LA ZONA QUE ESTAN EN LA FRONTERA PERO NO EN EL CONTORNO DEL SOLIDO COMPLETO (LOS DE LAS ARISTAS DE LAS CARAS DEL CONTORNO COMPLETO PERO QUE ESTAN EN CONTACTO CON OTRA CARA DE LA FRONTERA PERO QUE NO ES DEL CONTORNO DEL SOLIDO COMPLETO, NO LOS CONSIDERAMOS)
    gdl_frontera_no_contorno_solido_zona=[]

    for cara in caras_contorno_zona
        if findall(in(cara),caras_contorno_solido)==[]
            for nodo in caras_nodos_metam_exacto[cara]
                for gdl in 6*nodo-5:6*nodo
                    if findall(in(gdl), gdl_contorno_solido_completo_no_fronteras)==[] && findall(in(gdl), gdl_caras_repetidas)==[]
                        append!(gdl_frontera_no_contorno_solido_zona, gdl)
                    else
                    end
                end
            end
        else
        end
    end
    unique!(sort!(gdl_frontera_no_contorno_solido_zona))
    println(gdl_caras_repetidas);println()
    println(gdl_contorno_solido_completo_no_fronteras);println()
    println(gdl_frontera_no_contorno_solido_zona);println()



    ### OBTENER TODOS LOS GDL DE LA ZONA
    gdl_total_zona=[]

    for elemento_solido in elementos_alrededor
        for nodo in coords_metam_hexaedros_exacto[elemento_solido][:,4]
            append!(gdl_total_zona, collect(6*nodo-5:6*nodo))
        end
    end


    ### OBTENER LOS GDL FIJOS DE LA ZONA: LOS GDL RESTRINGIDOS DEL SOLIDO COMPLETO + TODOS LOS DE LA FRONTERA DE LA ZONA QUE NO SEAN DEL CONTORNO DEL SOLIDO COMPLETO
    gdl_fijos_zona=[]

    indices= findall(in(gdl_total_zona),gdl_restring_globales)
    gdl_fijos_zona=gdl_restring_globales[indices]

    append!(gdl_fijos_zona, gdl_frontera_no_contorno_solido_zona)
    unique!(sort!(gdl_fijos_zona))


    ### OBTENER TODOS LOS GDL DE LA ZONA NO FIJOS: TODOS LOS DE LAS CARAS DE DENTRO QUE NO ESTEN EN RESTRINGIDOS GLOBALES MÁS LOS DEL CONTORNO DEL SOLIDO REAL QUE ESTEN EN LAS ARISTAS ENTRE DOS ELEMENTOS
    gdl_no_fijos_zona=[]

    # por considerar los nodos de las caras, aqui ya se tiene en cuenta a los nodos de las aristas entre distintos hexaedros
    for cara in caras_repetidas_zona
        for nodo in caras_nodos_metam_exacto[cara]
            for gdl in 6*nodo-5:6*nodo
                if findall(in(gdl), gdl_fijos_zona)==[]
                    push!(gdl_no_fijos_zona, gdl)
                else
                end
            end
        end
    end
    unique!(sort!(gdl_no_fijos_zona))

    println(gdl_no_fijos_zona);println();println(gdl_fijos_zona)
    return elementos_alrededor, gdl_total_zona, gdl_no_fijos_zona, gdl_fijos_zona


end #function


function clasificacion_gdl_zona(elemento_solido, X,Y,Z,N_3D_x,N_3D_y,N_3D_z, coords_metam, conectiv_metam, N_capas_alrededor, hexaedros_alrededor_1_capa, caras_hexaedros, caras_contorno_solido, caras_nodos_metam_exacto, conectiv_metam_hexaedros_exacto, coords_metam_hexaedros_exacto, coords_metam_fronteras_exactas_hexaedros, gdl_restring_globales, gdl_contorno_solido_completo_no_fronteras, tensiones_pandeo_hexaedros)
    

    ### OBTENER LOS ELEMENTOS QUE ESTAN ALREDEDOR DEL QUE ESTUDIAMOS (TENIENDO EN CUENTA EL NUMERO DE CAPAS): INCLUIMOS AL PROPIO ELEMENTO
    elementos_alrededor=[elemento_solido]
    if N_capas_alrededor=="todas"
        N_capas_alrededor=maximum([N_3D_x, N_3D_y, N_3D_z])
    else
    end

    for capa in 1:N_capas_alrededor
        elementos_alrededor_nuevo=[]

        for elemento in elementos_alrededor
            #println(elemento)
            #println(hexaedros_alrededor_1_capa[elemento])
            append!(elementos_alrededor_nuevo, hexaedros_alrededor_1_capa[elemento])
            #readline()
        end

        append!(elementos_alrededor, elementos_alrededor_nuevo)
        unique!(sort!(elementos_alrededor))
    end


    ### OBTENER LA MATRIZ DE COORDENADAS DE LA ZONA COMPLETA
    nodos_zona=[]
    for elemento_solido in elementos_alrededor
        append!(nodos_zona, coords_metam_hexaedros_exacto[elemento_solido][:,4])
    end
    unique!(sort!(nodos_zona))
    coords_zona=coords_metam[nodos_zona,:]
    
    ### OBTENER LA MATRIZ DE CONECTIVIDAD DE LA ZONA COMPLETA
    barras_zona=[]; tensiones_pandeo_zona=Dict()
    for elemento_solido in elementos_alrededor
        for i in 1:length(conectiv_metam_hexaedros_exacto[elemento_solido][:,3])
            barra=conectiv_metam_hexaedros_exacto[elemento_solido][i,3]
            if findall(in(barra),barras_zona)==[]
                append!(barras_zona, barra)
                tensiones_pandeo_zona[barra]=tensiones_pandeo_hexaedros[elemento_solido][barra]
            else
            end
        end
    end
    unique!(sort!(barras_zona))
    conectiv_zona=conectiv_metam[barras_zona,:]    



    ### OBTENER LOS GDL ASOCIADOS AL CONTORNO DE ESTA ZONA
    caras_totales_zona=[]; caras_repetidas_zona=[]; caras_contorno_zona=[]

    # OBTENER LAS CARAS QUE SON DEL CONTORNO DE LA ZONA (SERAN LAS QUE NO SE REPITAN)
    for elemento in elementos_alrededor
        
        for cara in caras_hexaedros[elemento,:]
            if findall(in(cara), caras_totales_zona)==[]
                append!(caras_totales_zona, cara)
            else
                append!(caras_repetidas_zona, cara)
            end
        end

        indices_caras_no_repetidas=findall(!in(caras_repetidas_zona), caras_totales_zona)
        caras_contorno_zona=caras_totales_zona[indices_caras_no_repetidas]
    end

    # OBTENER LOS NODOS Y GDL DE ESA FRONTERA
    gdl_contorno_zona=[]

    for cara in caras_contorno_zona
        for nodo in caras_nodos_metam_exacto[cara]
            append!(gdl_contorno_zona, collect(6*nodo-5:6*nodo))
        end
    end
    unique!(sort!(gdl_contorno_zona))

    # OBTENER LOS GDL CORRESPONDIENTES A LAS CARAS QUE SE REPITEN
    gdl_caras_repetidas=[]
    for cara in caras_repetidas_zona
        for nodo in caras_nodos_metam_exacto[cara]
            append!(gdl_caras_repetidas, collect(6*nodo-5:6*nodo))
        end
    end


    ### OBTENER TODOS LOS GDL DE LA ZONA QUE ESTAN EN LA FRONTERA PERO NO EN EL CONTORNO DEL SOLIDO COMPLETO (LOS DE LAS ARISTAS DE LAS CARAS DEL CONTORNO COMPLETO PERO QUE ESTAN EN CONTACTO CON OTRA CARA DE LA FRONTERA PERO QUE NO ES DEL CONTORNO DEL SOLIDO COMPLETO, NO LOS CONSIDERAMOS)
        # este "if" es el que diferencia entre coger los desplazamientos del solido para el contorno tambien o dejarlo libre
    gdl_frontera_no_contorno_solido_zona=[]

    for cara in caras_contorno_zona
        if findall(in(cara),caras_contorno_solido)==[]
            for nodo in caras_nodos_metam_exacto[cara]
                for gdl in 6*nodo-5:6*nodo
                    if  findall(in(gdl), gdl_caras_repetidas)==[] && findall(in(gdl), gdl_contorno_solido_completo_no_fronteras)==[]
                        append!(gdl_frontera_no_contorno_solido_zona, gdl)
                    else
                    end
                end
            end
        else
        end
    end
    unique!(sort!(gdl_frontera_no_contorno_solido_zona))
    #println(gdl_caras_repetidas);println()
    #println(gdl_contorno_solido_completo_no_fronteras);println()
    #println(gdl_frontera_no_contorno_solido_zona);println()



    ### OBTENER TODOS LOS GDL DE LA ZONA
    gdl_total_zona=[]

    for elemento_solido in elementos_alrededor
        for nodo in coords_metam_hexaedros_exacto[elemento_solido][:,4]
            append!(gdl_total_zona, collect(6*nodo-5:6*nodo))
        end
    end
    unique!(sort!(gdl_total_zona))


    ### OBTENER LOS GDL FIJOS DE LA ZONA: LOS GDL RESTRINGIDOS DEL SOLIDO COMPLETO + TODOS LOS DE LA FRONTERA DE LA ZONA QUE NO SEAN DEL CONTORNO DEL SOLIDO COMPLETO
    gdl_fijos_zona=[]

    indices= findall(in(gdl_total_zona),gdl_restring_globales)
    gdl_fijos_zona=gdl_restring_globales[indices]

    append!(gdl_fijos_zona, gdl_frontera_no_contorno_solido_zona)
    unique!(sort!(gdl_fijos_zona))


    #println(gdl_fijos_zona)
    return elementos_alrededor, gdl_total_zona, gdl_fijos_zona, coords_zona, conectiv_zona, tensiones_pandeo_zona


end #function


function clasificacion_gdl_zona_mallado_fino_grueso(elemento_solido, coords_metam, conectiv_metam, conectiv_metam_hexaedros_exacto, coords_metam_hexaedros_exacto, coords_metam_fronteras_exactas_hexaedros, gdl_restring_globales, gdl_contorno_solido_completo_no_fronteras)
    

    ### OBTENER TODOS LOS GDL DE LA ZONA QUE ESTAN EN LA FRONTERA PERO NO EN EL CONTORNO DEL SOLIDO COMPLETO (LOS DE LAS ARISTAS DE LAS CARAS DEL CONTORNO COMPLETO PERO QUE ESTAN EN CONTACTO CON OTRA CARA DE LA FRONTERA PERO QUE NO ES DEL CONTORNO DEL SOLIDO COMPLETO, NO LOS CONSIDERAMOS)
    # este "if" es el que diferencia entre coger los desplazamientos del solido para el contorno tambien o dejarlo libre
    gdl_frontera_no_contorno_solido_zona=[]

    for nodo in coords_metam_fronteras_exactas_hexaedros[elemento_solido][:,4]
        for gdl in 6*nodo-5:6*nodo           
            if  findall(in(gdl), gdl_contorno_solido_completo_no_fronteras)==[]
                append!(gdl_frontera_no_contorno_solido_zona, gdl)
            else
            end
        end
    end

    unique!(sort!(gdl_frontera_no_contorno_solido_zona))


    ### OBTENER TODOS LOS GDL DE LA ZONA
    gdl_total_zona=[]

    for nodo in coords_metam_hexaedros_exacto[elemento_solido][:,4]
        append!(gdl_total_zona, collect(6*nodo-5:6*nodo))
    end

    unique!(sort!(gdl_total_zona))


    ### OBTENER LOS GDL FIJOS DE LA ZONA: LOS GDL RESTRINGIDOS DEL SOLIDO COMPLETO + TODOS LOS DE LA FRONTERA DE LA ZONA QUE NO SEAN DEL CONTORNO DEL SOLIDO COMPLETO
    gdl_fijos_zona=[]

    indices= findall(in(gdl_total_zona),gdl_restring_globales)
    gdl_fijos_zona=gdl_restring_globales[indices]

    append!(gdl_fijos_zona, gdl_frontera_no_contorno_solido_zona)
    unique!(sort!(gdl_fijos_zona))


    #println(gdl_fijos_zona)
    return gdl_total_zona, gdl_fijos_zona


end #function


function clasificacion_gdl_zona_mallado_fino_grueso_gradiente(elemento_solido, coords_metam, conectiv_metam, conectiv_metam_hexaedros_exacto, coords_metam_hexaedros_exacto, coords_metam_fronteras_exactas_hexaedros, gdl_restring_globales, gdl_contorno_solido_completo_no_fronteras)
    

    ### OBTENER TODOS LOS GDL DE LA ZONA QUE ESTAN EN LA FRONTERA PERO NO EN EL CONTORNO DEL SOLIDO COMPLETO (LOS DE LAS ARISTAS DE LAS CARAS DEL CONTORNO COMPLETO PERO QUE ESTAN EN CONTACTO CON OTRA CARA DE LA FRONTERA PERO QUE NO ES DEL CONTORNO DEL SOLIDO COMPLETO, NO LOS CONSIDERAMOS)
    # este "if" es el que diferencia entre coger los desplazamientos del solido para el contorno tambien o dejarlo libre
    gdl_frontera_no_contorno_solido_zona=[]

    for nodo in coords_metam_fronteras_exactas_hexaedros[elemento_solido][:,4]
        for gdl in 6*nodo-5:6*nodo           
            if  findall(in(gdl), gdl_contorno_solido_completo_no_fronteras)==[]
                append!(gdl_frontera_no_contorno_solido_zona, gdl)
            else
            end
        end
    end

    unique!(sort!(gdl_frontera_no_contorno_solido_zona))


    ### OBTENER TODOS LOS GDL DE LA ZONA
    gdl_total_zona=[]

    for nodo in coords_metam_hexaedros_exacto[elemento_solido][:,4]
        append!(gdl_total_zona, collect(6*nodo-5:6*nodo))
    end

    unique!(sort!(gdl_total_zona))


    ### OBTENER LOS GDL FIJOS DE LA ZONA: LOS GDL RESTRINGIDOS DEL SOLIDO COMPLETO + TODOS LOS DE LA FRONTERA DE LA ZONA QUE NO SEAN DEL CONTORNO DEL SOLIDO COMPLETO
    gdl_fijos_zona=[]

    indices= findall(in(gdl_total_zona),gdl_restring_globales)
    gdl_fijos_zona=gdl_restring_globales[indices]


    ### OBTENER LOS GDL NO FIJOS DE LA ZONA: TODOS LOS DE LA FRONTERA DE LA ZONA QUE NO SEAN DEL CONTORNO DEL SOLIDO COMPLETO
    gdl_no_fijos_zona=[]

    for gdl in gdl_frontera_no_contorno_solido_zona
        if findall(in(gdl),gdl_fijos_zona)==[]
            push!(gdl_no_fijos_zona, gdl)
        else
        end
    end
    unique!(sort!(gdl_no_fijos_zona))


    #println(gdl_fijos_zona)
    return gdl_total_zona, gdl_fijos_zona, gdl_no_fijos_zona


end #function


function resolucion_iterativa_acoplamiento_zonas_mallado_fino_grueso(hiperparametros ,X,Y,Z, N_3D_x_grueso, N_3D_y_grueso, N_3D_z_grueso, hexaedros_analizados_metam, hexaedros_mallado_fino_en_grueso, nodos_metam_en_hexaedros_finos, u_solido_acumulado, f_solido_acumulado, f_acumulado_hexaedros, energias_hexaedros_acumulado, epsilon_voigt_acumulado, f_acumulada_elementos_memoria, energias_acumuladas_elementos_memoria, eps_x_ptos_int_acumulado_memoria, eps_xy_ptos_int_acumulado_memoria, eps_xz_ptos_int_acumulado_memoria, coords_metam, conectividad_metam, coords_3D_fino, conectividad_3D_fino, conectividad_3D_grueso, volumenes_hexaedros_grueso, Youngs_elementos_finos_3D, Youngs_tangentes_hexaedros_gruesos, Youngs_tangentes_memoria_metam, tensiones_pandeo_elementos, punto_restriccion, gdl_restriccion, punto_desplazamiento, inc_valor_desplazamiento, coords_metam_hexaedros_exacto, conectiv_metam_hexaedros_exacto, porcentajes_barras_dentro_hexaedros, coords_metam_fronteras_exactas_hexaedros, gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam, gdl_total_zona_hexaedros, gdl_fijos_zona_hexaedros, gdl_total_hexaedro, gdl_fijos_hexaedro, gdl_restring_totales_hexaedro_local, gdl_contorno_solido_completo_no_fronteras, invariante_def_memoria_metam, iter_desplazamiento, radio)
    
    Young_solido_no_daño=curva_tension_deformacion_conjunto(10)

    # INICIALIZAR PARA EL ANALISIS DEL SOLIDO
    N_elementos_solidos_grueso=length(conectividad_3D_grueso[:,1])
    inc_u_solido=0; inc_f_solido=0; epsilon_voigt_nuevo=0; inc_energia_interna_total_solido=0
    young_material=Young_material(); SSR_antiguo=0
    hexaedros_analizados_metam_nuevo=copy(hexaedros_analizados_metam); # todos los elementos analizados en este desplazamiento
    inc_energia_hexaedros_total=0;
    energias_hexaedros_acumulado_nuevo=[]
    inc_deform_equiv_hexaed_aparente=zeros(N_elementos_solidos_grueso); inc_deform_equiv_solido_grueso=zeros(N_elementos_solidos_grueso); deform_equiv_solido_grueso=zeros(N_elementos_solidos_grueso) 
    inc_energias_solidos_gruesos=zeros(N_elementos_solidos_grueso); energias_solidos_gruesos=zeros(N_elementos_solidos_grueso)
    u_solido_acumulado_nuevo=[]; f_solido_acumulado_nuevo=[]; fuerzas_hexaedros_acumulado_nuevo=[]
    inc_energias_hexaedros_grueso_aparente=zeros(N_elementos_solidos_grueso)
    Youngs_tangentes_hexaedros_finos=copy(Youngs_elementos_finos_3D)
    
    
    # INICIALIZAR PARA EL ANALISIS DEL METAMATERIAL
    inc_u_metam=Dict(); inc_f_metam=Dict(); inc_energia_zona_metam=Dict()
    inc_energias_elementos=Dict(); inc_energia_elemento_solido_grueso_metam=Dict(); energia_elemento_solido_grueso_metam=Dict()
    Youngs_tangentes_memoria_metam_nuevo=Dict()
    f_acumulada_elementos_nuevo=Dict(); energias_acumuladas_elementos_nuevo=Dict()
    eps_x_ptos_int_nuevo=Dict(); eps_xy_ptos_int_nuevo=Dict(); eps_xz_ptos_int_nuevo=Dict()
    invariante_def_memoria_metam_nuevo=Dict()

    
    
    # HIPERPARAMETROS
    iter_max=hiperparametros[2];  tolerancia=hiperparametros[3]; tolerancia_2=hiperparametros[4]; diferencial=hiperparametros[5]
    eps_limite=hiperparametros[1]; modelo=hiperparametros[6]

    ### si el exponente es cero, entonces, se les pone a todos los hexaedros finos el mismo Young que al hexaedro grueso
    exponente_reparto_young=hiperparametros[7]

    SSR_vector=[]  # Vector error entre iteraciones


    # ITERAR PARA CONVERGER EN LA SOLUCION DEL MULTIESCALA
    for iter in 1:iter_max

        ######################### SOLIDO CONTINUO (MALLADO FINO) ######################################
        println("Iteracion desplazamiento ITER: ",iter ,"/", iter_max)
        println("### Inicio analisis del solido continuo (mallado fino)")

        ### TRADUCIR LAS CONDICIONES DE CONTORNO DEL SOLIDO A GRADOS DE LIBERTAD DEL SOLIDO
        gdl_restring_globales_solido, inc_valor_restricciones_solido=condic_contorno_solido(punto_restriccion, gdl_restriccion, punto_desplazamiento, inc_valor_desplazamiento, coords_3D_fino)

       
        ### RESOLVER LOS INCREMENTOS DE DESPLAZAMIENTOS DEL SOLIDO        
        inc_u_solido, inc_f_solido, inc_energias_hexaedros, deform_equiv_solido, inc_deform_equiv_solido, epsilon_voigt_nuevo = resolucion_solido_3D(coords_3D_fino, conectividad_3D_fino, Youngs_tangentes_hexaedros_finos, epsilon_voigt_acumulado, gdl_restring_globales_solido, inc_valor_restricciones_solido)
        
        ### INCREMENTO DE ENERGIA TOTAL DEL SOLIDO
        inc_energia_interna_total_solido=sum(inc_energias_hexaedros)        

        ### DESPLAZAMIENTOS Y FUERZAS TOTALES A PARTIR DE LOS INCREMENTOS 
        u_solido_acumulado_nuevo = u_solido_acumulado + inc_u_solido
        f_solido_acumulado_nuevo = f_solido_acumulado + inc_f_solido

        ### ENERGIA TOTAL ACUMULADA EN CADA ELEMENTO Y VECTOR DE FUERZAS EN CADA ELEMENTO
        energias_hexaedros_acumulado_nuevo, fuerzas_hexaedros_acumulado_nuevo = energias_y_fuerzas_acumuladas_elementos(energias_hexaedros_acumulado, inc_energias_hexaedros, f_acumulado_hexaedros, inc_u_solido, coords_3D_fino, conectividad_3D_fino, Youngs_tangentes_hexaedros_finos)
        
        #$$$
        #println();println("ITER: ",iter);
        println("deformacion equiv. max solido homogeneo= ",maximum(deform_equiv_solido));println("inc. energia solido homogeneo= ", sum(inc_energias_hexaedros))#;println("Youngs= ", Youngs_tangentes_hexaedros_finos);println()
        #println();println("Youngs_hexaed");println(Youngs_tangentes_hexaedros_finos);println(u_solido);println(fuerzas_nodos_solido)
        #readline()


        ############################ ENERGIAS Y DEFORMACIONES DE LOS SOLIDOS GRUESOS EN FUNCION DE LOS FINOS ##############################
        println("### inicio analisis del solido continuo con mallado grueso (hexaedro grueso como suma de hexaedros finos)")

        for elemento_solido_grueso in 1:N_elementos_solidos_grueso  

            inc_energias_solidos_gruesos[elemento_solido_grueso]=0
            energias_solidos_gruesos[elemento_solido_grueso]=0
            deform_equiv_solido_grueso[elemento_solido_grueso]=0; suma_2=0
            inc_deform_equiv_solido_grueso[elemento_solido_grueso]=0; suma_1=0

            N_hexaedros_finos_en_grueso=length(hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:])

            for elemento_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                inc_energias_solidos_gruesos[elemento_solido_grueso]=inc_energias_solidos_gruesos[elemento_solido_grueso] + inc_energias_hexaedros[elemento_fino]
                energias_solidos_gruesos[elemento_solido_grueso]=energias_solidos_gruesos[elemento_solido_grueso] + energias_hexaedros_acumulado_nuevo[elemento_fino]
                suma_1=suma_1 + inc_deform_equiv_solido[elemento_fino]^2/N_hexaedros_finos_en_grueso
                suma_2=suma_2 + deform_equiv_solido[elemento_fino]^2/N_hexaedros_finos_en_grueso
            end
            inc_deform_equiv_solid=sqrt(suma_1)#; println("1 ",sqrt(suma_1))
            deform_equiv_solido_grueso[elemento_solido_grueso]=sqrt(suma_2)
            inc_deform_equiv_solido_grueso[elemento_solido_grueso]=sqrt(inc_energias_solidos_gruesos[elemento_solido_grueso]/(1/2*Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]*volumenes_hexaedros_grueso[elemento_solido_grueso]))
            #println("2 ", inc_deform_equiv_solido_grueso[elemento_solido_grueso])
        end
        #readline()

        ################################ ACOPLAMIENTO MULTIESCALA ####################################

        ### PASAR LAS RESTRICCIONES DEL SOLIDO A RESTRICCIONES DEL METAMATERIAL (CON 6 GDL EN VEZ DE 3)
        println("### inicio acoplamiento solido-metamaterial: relacion entre GDLs y puntos")

        if gdl_restring_globales_metam==[]
            N_restricciones = length(gdl_restriccion[:,1])
            N_desplazamientos = length(inc_valor_desplazamiento[:,1])

            punto_restriccion_metam = punto_restriccion
            gdl_restriccion_metam = [gdl_restriccion zeros(N_restricciones,3)]
            punto_desplazamiento_metam = punto_desplazamiento
            inc_valor_desplazamiento_metam = [inc_valor_desplazamiento zeros(N_desplazamientos,3)]

            # GRADOS DE LIBERTAD CON ALGUNA RESTRICCION Y EL VALOR DE ESA RESTRICCION (PUEDEN SER DEL CONTORNO O NO) del solido completo
            gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam=condic_contorno(punto_restriccion_metam, gdl_restriccion_metam ,punto_desplazamiento_metam ,inc_valor_desplazamiento_metam ,coords_metam)
            #println(gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam)

            
            # OBTENER LOS GDL DEL CONTORNO COMPLETO QUE NO ESTEN EN LAS FRONTERAS EN CONTACTO CON OTRAS ZONAS
            Nodos_metam_contorno_solido_no_fronteras, gdl_contorno_solido_completo_no_fronteras , Nodos_meetam_contorno_solido_fronteras, gdls_contorno_solido_completo_fronteras=gdl_metamaterial_contorno_solido_no_fronteras(X,Y,Z,N_3D_x_grueso, N_3D_y_grueso, N_3D_z_grueso, coords_metam)

        else
        end


        ### SI EN ALGUN ELEMENTO SOLIDO DEL MALLADO FINO SE SOBREPASA LA TENSION LIMITE (ARBITRARIA) SE HACE EL ANALISIS MULTIESCALA DEL ELEMENTO DEL MALLADO GRUESO QUE CONTENGA A ESE ELEMENTO DEL MALLADO FINO
        println("### Analisis multiescala de elementos solidos que han superado el umbral")

        for elemento_solido_grueso in 1:N_elementos_solidos_grueso  
    
            println();println("ELEMENTO hexaedrico grueso ",elemento_solido_grueso, " Y SU ZONA")#;println("deform_equiv_tet= ", deform_equiv_solido[elemento_solido_grueso]);println()
            
            # COMPROBAR SI EN ALGUN ELEMENTO FINO DE LOS QUE COMPONEN AL GRUESO SE SUPERA LA DEFORMACION LIMITE
            valor=0
            for elemento_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                if deform_equiv_solido[elemento_fino]<= eps_limite
                else
                    valor=1
                    break
                end
            end

            # si no se supera el limite, el young del elemento será el de referencia del solido compacto
            if valor==0 
                println("El hexaedro no supera el umbral de deformacion")
            
            # si se supera un valor umbral, se calcula el valor del young equivalente de esa zona con el modelo multiescala
            else   
                println("El hexaedro supera el umbral de deformacion --> se estudia con metamaterial")
                push!(hexaedros_analizados_metam_nuevo, elemento_solido_grueso)

                ### GENERAR LOS GDL DE LA ZONA
                # SI NO SE HA ESTUDIADO TODAVIA EL HEXAEDRO (Y SU ZONA CORRESPONDIENTE ALREDEDOR, SE GENERAN LOS GDL CORRESPONDIENTES)
                # COMO NO RELAJAMOS, ENTONCES LOS GDL FIJOS SON IGUALES A LOS RESTRINGIDOS
                if haskey(gdl_total_zona_hexaedros, elemento_solido_grueso)== false
                    gdl_total_zona_hexaedros[elemento_solido_grueso], gdl_fijos_zona_hexaedros[elemento_solido_grueso] =clasificacion_gdl_zona_mallado_fino_grueso(elemento_solido_grueso, coords_metam, conectividad_metam, conectiv_metam_hexaedros_exacto, coords_metam_hexaedros_exacto, coords_metam_fronteras_exactas_hexaedros, gdl_restring_globales_metam, gdl_contorno_solido_completo_no_fronteras)

                    # GDL CON RESTRICCIONES FIJAS DE LA ZONA (SON TODOS LOS GDL RESTRINGIDOS)
                    gdl_restring_totales_hexaedro_local[elemento_solido_grueso]=indexin(gdl_fijos_zona_hexaedros[elemento_solido_grueso], gdl_total_zona_hexaedros[elemento_solido_grueso])

                else
                end
                

                # SI TODAVIA NO SE HA ESTUDIADO ESTA ZONA, SE PONE LA CONDICION DE CONTORNO DEL DESPLAZAMIENTO TOTAL DEL SOLIDO
                if findall(in(elemento_solido_grueso), hexaedros_analizados_metam)==[] #si no ha sido estudiada esta zona para desplazamientos iteraciones anteriores
                    println("nuevo elemento estudiado con multiescala en este desplazamiento")
                    inc_u_solido_condic_contorno_multiescala=u_solido_acumulado_nuevo
                    inc_valor_restricciones_fijos_global_metam_aparente=copy(inc_valor_restricciones_fijos_global_metam)*iter_desplazamiento

                # SI SE HA ESTUDIADO YA ESTA ZONA, LE IMPONEMOS EL INCREMENTO DE DESPLAZAMIENTO DEL SOLIDO
                else
                    inc_u_solido_condic_contorno_multiescala=inc_u_solido
                    inc_valor_restricciones_fijos_global_metam_aparente=copy(inc_valor_restricciones_fijos_global_metam)
                end


                ### OBTENER LAS CONDICIONES DE CONTORNO PARA LOS NODOS RESTRINGIDOS FIJOS UTILIZANDO TANTO LOS DESPLAZAMIENTOS DEL SOLIDO COMO LAS RESTRICCIONES GLOBALES
                N_gdl_restricciones_fijas_zona=length(gdl_fijos_zona_hexaedros[elemento_solido_grueso])
                inc_valor_restricciones_fijas_zona=zeros(N_gdl_restricciones_fijas_zona)

                for i in 1:N_gdl_restricciones_fijas_zona
                    gdl=gdl_fijos_zona_hexaedros[elemento_solido_grueso][i]
                    #println(gdl)

                    # Si el gdl esta en las restricciones globales
                    if (indice=findall(in(gdl),gdl_restring_globales_metam))!=[]
                        inc_valor_restricciones_fijas_zona[i]=inc_valor_restricciones_fijos_global_metam_aparente[indice[1]]

                    # Si el gdl no esta en las restricciones globales, se le restringe a partir de los desplazamientos del solido
                    else
                        # DETERMINAR A QUE NODO LE CORRESPONDE CADA GDL Y QUE VALOR DE DESPLAZAMIENTO LE TOCA
                        nodo=Int64(ceil(gdl/6))
                        #elementos_nodo=0; inc_desplaz=0

                        for elemento_solido_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                            indice=findall(in(nodo), nodos_metam_en_hexaedros_finos[elemento_solido_fino])
                            #elementos_nodo=elementos_nodo+1                       
                            
                            # SI ESTÁ EN EL ELEMENTO SOLIDO DE LA ZONA QUE ESTAMOS MIRANDO AHORA
                            if indice!=[]
                                indice_nodo=findall(in(nodo),coords_metam_fronteras_exactas_hexaedros[elemento_solido_grueso][:,4])
                                coords_punto=coords_metam_fronteras_exactas_hexaedros[elemento_solido_grueso][indice_nodo[1],1:3]'
                                
                                resto=gdl%6 # como imponemos desplazamientos y giros, el resto irá entra 0 y 5

                                # SI EL GDL CORRESPONDE A UN DESPLAZAMIENTO EN X
                                if resto == 1
                                    inc_desplaz_nodo=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                                    
                                    #inc_valor_restricciones_fijas_zona[i]=inc_desplaz_nodo
                                    #break

                                # SI EL GDL CORRESPONDE A UN DESPLAZAMIENTO EN Y    
                                elseif resto == 2
                                    inc_desplaz_nodo=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                                            
                                    #inc_valor_restricciones_fijas_zona[i]=inc_desplaz_nodo
                                    #break

                                # SI EL GDL CORRESPONDE A UN DESPLAZAMIENTO EN Z    
                                elseif resto == 3
                                    inc_desplaz_nodo=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                                            
                                    #inc_valor_restricciones_fijas_zona[i]=inc_desplaz_nodo
                                    #break

                                # SI EL GDL CORRESPONDE A UN GIRO EN X (AROXIMAMOS EL GIRO)   
                                elseif resto == 4                        
                                    inc_desplaz_nodo_y=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                                    inc_desplaz_nodo_diferencial_y=desplazamientos_elemento(coords_punto+[0 0 diferencial], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                                    inc_desplaz_nodo_1=-(inc_desplaz_nodo_diferencial_y-inc_desplaz_nodo_y)/ diferencial

                                    inc_desplaz_nodo_z=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                                    inc_desplaz_nodo_diferencial_z=desplazamientos_elemento(coords_punto+[0 diferencial 0], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                                    inc_desplaz_nodo_2=(inc_desplaz_nodo_diferencial_z-inc_desplaz_nodo_z)/ diferencial
             
                                    inc_desplaz_nodo=(inc_desplaz_nodo_1+inc_desplaz_nodo_2)/2
                                    #inc_desplaz=inc_desplaz+(inc_desplaz_nodo_1+inc_desplaz_nodo_2)
                                    #inc_desplaz_nodo=inc_desplaz/(2*elementos_nodo)
                                    
                                # SI EL GDL CORRESPONDE A UN GIRO EN Y (AROXIMAMOS EL GIRO)   
                                elseif resto == 5                        
                                    inc_desplaz_nodo_z=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                                    inc_desplaz_nodo_diferencial_z=desplazamientos_elemento(coords_punto+[diferencial 0 0], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                                    inc_desplaz_nodo_1=-(inc_desplaz_nodo_diferencial_z-inc_desplaz_nodo_z)/ diferencial
                                    
                                    inc_desplaz_nodo_x=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                                    inc_desplaz_nodo_diferencial_x=desplazamientos_elemento(coords_punto+[0 0 diferencial], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                                    inc_desplaz_nodo_2=(inc_desplaz_nodo_diferencial_x-inc_desplaz_nodo_x)/ diferencial
                                    
                                    inc_desplaz_nodo=(inc_desplaz_nodo_1+inc_desplaz_nodo_2)/2
                                    #inc_desplaz=inc_desplaz+(inc_desplaz_nodo_1+inc_desplaz_nodo_2)
                                    #inc_desplaz_nodo=inc_desplaz/(2*elementos_nodo)

                                # SI EL GDL CORRESPONDE A UN GIRO EN Z (AROXIMAMOS EL GIRO)   
                                elseif resto == 0                     
                                    inc_desplaz_nodo_y=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                                    inc_desplaz_nodo_diferencial_y=desplazamientos_elemento(coords_punto+[diferencial 0 0], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                                    inc_desplaz_nodo_1=(inc_desplaz_nodo_diferencial_y-inc_desplaz_nodo_y)/ diferencial
                                    
                                    inc_desplaz_nodo_x=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                                    inc_desplaz_nodo_diferencial_x=desplazamientos_elemento(coords_punto+[0 diferencial 0 ], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                                    inc_desplaz_nodo_2=-(inc_desplaz_nodo_diferencial_x-inc_desplaz_nodo_x)/ diferencial

                                    inc_desplaz_nodo=(inc_desplaz_nodo_1+inc_desplaz_nodo_2)/2
                                    #inc_desplaz=inc_desplaz+(inc_desplaz_nodo_1+inc_desplaz_nodo_2)
                                    #inc_desplaz_nodo=inc_desplaz/(2*elementos_nodo)
                                end

                                inc_valor_restricciones_fijas_zona[i]=inc_desplaz_nodo
                                # #println(gdl," ",inc_desplaz_nodo," ",elemento_solido_fino) 
                                break
                            else
                            end
                        end
                    end
                end

                #println(inc_valor_restricciones_fijas_zona)

                ################################# RESOLUCION DE LA ZONA ##################################
                println("### Resolucion con metamaterial de la zona")

                ########### OBTENER LOS GDL QUE VAN A INTERVENIR Y LA CONECTIVIDAD DENTRO DEL HEXAEDRO
                # coordenadas de los nodos de la zona con metamaterial (un elemento hexaedrico)
                coords_metam_elemento=coords_metam_hexaedros_exacto[elemento_solido_grueso]

                # Conectividad de las barras de la zona
                conectividad_zona_metam=conectiv_metam_hexaedros_exacto[elemento_solido_grueso]
                #porcentajes_barras_zona=porcentajes_barras_dentro_hexaedros[elemento_solido_zona]
                N_elementos_metam=length(conectividad_zona_metam[:,1])
                #$$$$
                #println();println("nodos_tot_extend");println(coords_metam_elemento_extendido);println()


                ########### INICIALIZAR VARIABLES Y DECIDIR QUE CONDICIONES DE CONTORNO
                # INICIALIZAR LA MATRIZ DE RIGIDEZ DE LA ZONA DE METAMATERIAL QUE ESTAMOS ESTUDIANDO

                # SI NO HA SIDO ESTUDIADA ESTA ZONA PARA DESPLAZAMIENTOS ITERACIONES ANTERIORES
                if findall(in(elemento_solido_grueso), hexaedros_analizados_metam)==[] 
                    
                    # si no hay registros anteriores, inicializamos las matrices
                    inicializar_cero=Dict(); eps_x_ptos_int_acumulado=Dict(); eps_xy_ptos_int_acumulado=Dict(); eps_xz_ptos_int_acumulado=Dict(); f_acumulada_elementos_zona=Dict() ;energias_acumuladas_elementos_zona=Dict(); invariante_def_metam=Dict()
                    for i in 1:N_elementos_metam
                        element_metam=conectividad_zona_metam[i,3]
                        inicializar_cero[element_metam]=0 # un diccionario inicializado a cero para uso generico
                        eps_x_ptos_int_acumulado[element_metam]=[0,0]      
                        eps_xy_ptos_int_acumulado[element_metam]=[0,0]
                        eps_xz_ptos_int_acumulado[element_metam]=[0,0]
                        energias_acumuladas_elementos_zona[element_metam]=0 
                        f_acumulada_elementos_zona[element_metam] =zeros(12)
                        invariante_def_metam[element_metam]=0    
                    end
                    Youngs_tangentes_in=curva_tension_deformacion_tabla_Osgood(inicializar_cero,"deform_equiv",inicializar_cero)  #si no se tienen registros de modulos de young anteriores, se empieza por el valor base del material
                    Youngs_tangentes_memoria_metam[elemento_solido_grueso]=Youngs_tangentes_in
                    inc_deform_equiv_hexaed_aparente[elemento_solido_grueso]=deform_equiv_solido_grueso[elemento_solido_grueso]
                    inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]=energias_solidos_gruesos[elemento_solido_grueso]
                    
                # SI LA ZONA YA HA SIDO ESTUDIADA, COGEMOS LOS VALORES QUE YA EXISTEN Y ESTUDIAMOS INCREMENTOS DESDE ESE PUNTO
                else                       
                    #si ya hay un vector de Youngs calculado para el desplazamiento anterior, se coje ese para facilitar la convergencia (cambio menos brusco)          
                    Youngs_tangentes_in=Dict(); eps_x_ptos_int_acumulado=Dict(); eps_xy_ptos_int_acumulado=Dict(); eps_xz_ptos_int_acumulado=Dict(); f_acumulada_elementos_zona=Dict() ;energias_acumuladas_elementos_zona=Dict(); invariante_def_metam=Dict()
                    
                    Youngs_tangentes_in=copy(Youngs_tangentes_memoria_metam[elemento_solido_grueso])
                    eps_x_ptos_int_acumulado=copy(eps_x_ptos_int_acumulado_memoria[elemento_solido_grueso])      
                    eps_xy_ptos_int_acumulado=copy(eps_xy_ptos_int_acumulado_memoria[elemento_solido_grueso])   
                    eps_xz_ptos_int_acumulado=copy(eps_xz_ptos_int_acumulado_memoria[elemento_solido_grueso])
                    energias_acumuladas_elementos_zona =copy(energias_acumuladas_elementos_memoria[elemento_solido_grueso])
                    f_acumulada_elementos_zona =copy(f_acumulada_elementos_memoria[elemento_solido_grueso])
                    invariante_def_metam=copy(invariante_def_memoria_metam[elemento_solido_grueso])

                    inc_deform_equiv_hexaed_aparente[elemento_solido_grueso]=inc_deform_equiv_solido_grueso[elemento_solido_grueso]
                    inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]=inc_energias_solidos_gruesos[elemento_solido_grueso]
                end
                                    



                ########### RESOLVER EL PROBLEMA DEL HEXAEDRO CON METAMATERIAL (EL ZOOM EN ESE ELEMENTO DEL SOLIDO)                 

                # obtener los incrementos en desplazamientos, la energia de la zona extedida y los modulos de Young
                inc_u_metam[elemento_solido_grueso], inc_f_metam[elemento_solido_grueso], f_acumulada_elementos_nuevo[elemento_solido_grueso], inc_energia_zona_metam[elemento_solido_grueso], inc_energias_elementos[elemento_solido_grueso], energias_acumuladas_elementos_nuevo[elemento_solido_grueso], Youngs_tangentes_memoria_metam_nuevo[elemento_solido_grueso], eps_x_ptos_int_nuevo[elemento_solido_grueso], eps_xy_ptos_int_nuevo[elemento_solido_grueso], eps_xz_ptos_int_nuevo[elemento_solido_grueso], invariante_def_memoria_metam_nuevo[elemento_solido_grueso]=resolucion_iterativa_metamaterial(coords_metam_elemento, conectividad_zona_metam, gdl_restring_totales_hexaedro_local[elemento_solido_grueso], inc_valor_restricciones_fijas_zona , Youngs_tangentes_in, tensiones_pandeo_elementos[elemento_solido_grueso], eps_x_ptos_int_acumulado, eps_xy_ptos_int_acumulado, eps_xz_ptos_int_acumulado, energias_acumuladas_elementos_zona, f_acumulada_elementos_zona, invariante_def_metam, radio)
                # println("desplazamientos y fuerzas")
                # for i in 1:Int64(length(inc_f_metam[elemento_solido_grueso])/6)
                #     println(i," ",inc_u_metam[elemento_solido_grueso][6*i-5:6*i]," ",inc_f_metam[elemento_solido_grueso, elemento_solido_grueso_zona][6*i-5:6*i])
                # end
                # println()
                # println("energia")
                # for i in conectividad_zona_metam[:,3]
                #     println(i," ",inc_energias_elementos[elemento_solido_grueso][i], " ", energias_acumuladas_elementos_nuevo[elemento_solido_grueso][i], " ", Youngs_tangentes_memoria_metam_nuevo[elemento_solido_grueso][i])
                #     indice=findall(in(i),conectividad_zona_metam[:,3])
                #     nodos=conectividad_zona_metam[indice[1],1:2]
                #     println(nodos)
                #     indices=findall(in(nodos),coords_metam_hexaedros_exacto[elemento_solido_grueso][:,4])
                #     println("   ",inc_u_metam[elemento_solido_grueso][6*indices[1]-5:6*indices[1]])
                #     println("   ",inc_u_metam[elemento_solido_grueso][6*indices[2]-5:6*indices[2]])
                # end

                #readline()


                ######## UNA VEZ ANALIZADOS TODOS LOS ELEMENTOS DE LA ZONA, SOLO NOS IMPORTA LA ENERGIA DEL ELEMENTO SOLIDO GRUESO

                # sumamos la energia de las barras del elemento
                inc_energia_hexaedro_metam=0; energia_hexaedro_metam=0; inc_energia_daño_cero=0
                for barra in conectiv_metam_hexaedros_exacto[elemento_solido_grueso][:,3]
                    inc_energia_hexaedro_metam=inc_energia_hexaedro_metam + inc_energias_elementos[elemento_solido_grueso][barra]*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]
                    energia_hexaedro_metam=energia_hexaedro_metam + energias_acumuladas_elementos_nuevo[elemento_solido_grueso][barra]*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]
                    
                    inc_energia_daño_cero=inc_energia_daño_cero+inc_energias_elementos[elemento_solido_grueso][barra]*young_material/Youngs_tangentes_memoria_metam_nuevo[elemento_solido_grueso][barra]*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]
                    #$$$$$
                    # indice=findall(in(barra),conectividad_metam[:,3])
                    # nodo_1=conectividad_metam[indice,1]; nodo_2=conectividad_metam[indice,2]
                    # println(barra," ", nodo_1," ", nodo_2, " ", porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]," ",inc_energias_elementos[elemento_solido_grueso][barra]*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra])
                end
                daño_hexaedro_grueso=inc_energia_hexaedro_metam/inc_energia_daño_cero    # parametro de daño definido como la energia exacta del metamaterial frente al valor de energia en el rango lineal  

                inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]=inc_energia_hexaedro_metam
                energia_elemento_solido_grueso_metam[elemento_solido_grueso]=energia_hexaedro_metam
                
                println();println("hexaedro= ",elemento_solido_grueso," inc_energia_exacto_metam= ",inc_energia_hexaedro_metam)
                #readline()


                ########### OBTENER EL PARAMETRO DE DAÑO QUE HEMOS DEFINIDO PARA EL REPARTO DEL YOUNG A CADA ELEMENTO FINO

                if (exponente_reparto_young==0) && (modelo != "modelo_daño_hexaedros_finos")
                else
                    daño_hexaedros_finos=Dict()
                    
                    for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                        inc_energia_hexaedros_finos=0; inc_energia_daño_nulo=0; barras=[]; daño_hexaedros_finos[hexaedro_fino]=0

                        for nodo in nodos_metam_en_hexaedros_finos[hexaedro_fino]
                            indices=findall(in(nodo), conectividad_zona_metam[:,1:2])
                            for indice in indices
                                push!(barras, conectividad_zona_metam[indice[1],3])
                            end
                            unique!(sort!(barras))
                        end
                        
                        for barra in barras                    
                            inc_energia_hexaedros_finos=inc_energia_hexaedros_finos + inc_energias_elementos[elemento_solido_grueso][barra]#*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]
                            inc_energia_daño_nulo=inc_energia_daño_nulo + inc_energias_elementos[elemento_solido_grueso][barra]*young_material/Youngs_tangentes_memoria_metam_nuevo[elemento_solido_grueso][barra]#*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]
                        end

                        daño_hexaedros_finos[hexaedro_fino]=inc_energia_hexaedros_finos/inc_energia_daño_nulo
                        #println(daño_hexaedros_finos[hexaedro_fino])
                    end
                end



                ########### OBTENER EL MODULO DE YOUNG TANGENTE NUEVO DE ESTE HEXAEDRO (APROXIMACION) --> IGUALAMOS LAS ENERGIAS, Y SUPONIENDO LOS DESPLAZAMIENTOS BIEN CALCULADOS, OBTENEMOS EL NUEVO YOUNG
                # Lo que hacemos es como calcular el modulo de Young tangente del hexaedro pero a partir de igualar la energia

                ### DEPENDIENDO DEL MODELO UTILIZADO SE ACTUALIZAN LOS YOUNGS TANGENTES DE UNA FORMA U OTRA
                println("### actualizacion del young tangente del solido continuo a partir de los calculos del metamaterial")

                if modelo == "modelo_energia_hexaedros_gruesos"    # Este modelo actualiza a todos los hexaedros finos con el mismo valor que el del hexaedro grueso. Los actualiza ajustando el Young para que la energía se iguale con la de la zona calculada por el metamaterial
                    
                    # Young_tangente_nuevo_2=inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]/(1/2*inc_deform_equiv_hexaed_aparente[elemento_solido_grueso]^2*volumenes_hexaedros_grueso[elemento_solido_grueso]) # nueva estimacion del valor del Young del hexaedro
                    Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]=inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]/(inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]/Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso])
                    
                    ### si el exponente es cero, entonces, se les pone a todos los hexaedros finos el mismo Young que al hexaedro grueso
                    if exponente_reparto_young==0
                        for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                            Youngs_tangentes_hexaedros_finos[hexaedro_fino]=Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]         
                        end
                    else
                        for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                            Youngs_tangentes_hexaedros_finos[hexaedro_fino]=Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]*(daño_hexaedros_finos[hexaedro_fino]/daño_hexaedro_grueso)^exponente_reparto_young          
                        end
                    end

                    println("Young hexaedro grueso modelo energia: ",Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso])#, " ", sum(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])/length(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]]))
                    println("Young hexaedro grueso modelo daño: ",Young_solido_no_daño*daño_hexaedro_grueso)#, " ", sum(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])/length(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]]))
                    println("youngs hexaedros finos: ", Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])
                    ##Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]].= Young_tangente_nuevo# (Youngs_elementos_finos_3D[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,1]] + Young_tangente_nuevo)/2
                
                elseif modelo == "modelo_daño_hexaedros_finos"    # Este modelo actualiza a todos los hexaedros finos con valores distintos a partir de el del hexaedro grueso. Los actualiza ajustando el Young para que la energía se iguale con la de la zona calculada por el metamaterial
                 
                    for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                        Youngs_tangentes_hexaedros_finos[hexaedro_fino]=Young_solido_no_daño*daño_hexaedros_finos[hexaedro_fino]             
                    end
                                       
                    Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]=Young_solido_no_daño*daño_hexaedro_grueso
                    println("Young hexaedro grueso modelo energia: ",inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]/(inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]/Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]))
                    println("Young hexaedro grueso modelo daño: ",Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso], " ", sum(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])/length(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]]))
                    println("youngs hexaedros finos: ", Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])
                
                elseif modelo == "modelo_daño_hexaedros_gruesos"   # Este modelo actualiza a todos los hexaedros finos con el mismo valor que el del hexaedro grueso. Los actualiza ajustando el Young en funcion del parametro de daño definido antes

                    if exponente_reparto_young==0
                        for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                            Youngs_tangentes_hexaedros_finos[hexaedro_fino]=Young_solido_no_daño*daño_hexaedro_grueso
                        end
                    else
                        for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                            Youngs_tangentes_hexaedros_finos[hexaedro_fino]=Young_solido_no_daño*daño_hexaedro_grueso*(daño_hexaedros_finos[hexaedro_fino]/daño_hexaedro_grueso)^exponente_reparto_young   
                        end
                    end
                    
                    Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]=Young_solido_no_daño*daño_hexaedro_grueso
                    println("Young hexaedro grueso modelo energia: ",inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]/(inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]/Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]))
                    println("Young hexaedro grueso modelo daño: ",Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso], " ", sum(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])/length(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]]))
                    println("youngs hexaedros finos: ", Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])
                
                else
                    println("No se ha introducido modelo de daño")
                end
                    #readline()
                #println(1/2*inc_deform_equiv_hexaed_aparente[elemento_solido_grueso]^2*volumenes_hexaedros_grueso[elemento_solido_grueso]," ", inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]/Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,1]])
                #println(Young_tangente_nuevo," ",Young_tangente_nuevo_2," ",Youngs_elementos_finos_3D[elemento_solido_grueso])    
                
                #$$$$
                println("inc. energ. hexaedro solido= ",inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]);println("inc. energ. metamaterial= ",inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]);println()
                println("energ. acumulado hexaedro solido= ",energias_solidos_gruesos[elemento_solido_grueso]);println("energia acumulada metamaterial= ",energia_elemento_solido_grueso_metam[elemento_solido_grueso]);println()
                #readline()

            end
            
        end
        
        println("Young's tangentes todos hexaedros gruesos", Youngs_tangentes_hexaedros_gruesos)
        unique!(sort!(hexaedros_analizados_metam_nuevo))

        ### RESUMEN DATOS ENERGIAS
        println("ENERGIA INCREMENTADA HEXAEDROS METAMATERIAL VS SOLIDO CONTINUO: ")
        for hexaedro in hexaedros_analizados_metam_nuevo 

            # si el elemento se estudia por primera vez en este incremento
            if (findall(in(hexaedro), hexaedros_analizados_metam_nuevo)!=[]) &&  (findall(in(hexaedro), hexaedros_analizados_metam)!=[])
                println("ELEMENTO: ", hexaedro, " inc energia como continuo= ",inc_energias_solidos_gruesos[hexaedro], " inc energia como metam= ",inc_energia_elemento_solido_grueso_metam[hexaedro])
            else
                println("ELEMENTO: ", hexaedro, " energia como continuo= ",energias_solidos_gruesos[hexaedro], " energia como metam= ",inc_energia_elemento_solido_grueso_metam[hexaedro])
            end
        end
        println()
        println("INC ENERGIA TOTAL SOLIDO CONTINUO ",sum(inc_energias_hexaedros))
        println("ENERGIA TOTAL SOLIDO CONTINUO ",sum(energias_hexaedros_acumulado_nuevo))

        # SI SE ANALIZAN TODOS LOS HEXAEDROS
        if length(hexaedros_analizados_metam_nuevo)==N_elementos_solidos_grueso
            println("INC ENERGIA TOTAL METAMATERIAL ",sum(values(inc_energia_elemento_solido_grueso_metam)))
            println("ENERGIA TOTAL METAMATERIAL ",sum(values(energia_elemento_solido_grueso_metam)))
        else
        end
        
        #readline()


        ### CRITERIO DE SALIDA
        # si no se ha analizado ningun hexaedro se sale directamente
        if hexaedros_analizados_metam_nuevo==[]
            break
        # si se han analizado hexaedros, vemos si se cumple el criterio de salida
        else
            SSR=0; N_hexaedros_analizados=length(hexaedros_analizados_metam_nuevo)
            for elemento_solido_grueso in hexaedros_analizados_metam_nuevo #collect(keys(inc_energia_elemento_solido_zona))
                # el criterio de salida será con sum of square residuals
                #SSR = SSR + 1/N_hexaedros_analizados*(inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]-inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso])^2
                SSR=SSR+abs(inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]-inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso])/inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]
                
            end

            push!(SSR_vector, SSR)
            println();println("error SSR= ",SSR)
            if (SSR<tolerancia)#(abs(sqrt(SSR)-sqrt(SSR_antiguo))/sqrt(SSR) < tolerancia) || (abs(sqrt(SSR))/inc_energia_hexaedros_total<tolerancia_2)
                # Youngs_tangentes_hexaedros_finos=Youngs_tangentes_hexaedros_finos # este vector con los valores de los modulos de young lo usaremos para el siguiente valor de desplazamiento impuesto y asi facilitar la convergencia
                # println(Youngs_tangentes_hexaedros_finos)
                break
            else; 
                SSR_antiguo=SSR #;readline()
                if (iter==iter_max) && (abs(SSR) != minimum(abs.(SSR_vector)))
                    println("ATENCION: No ha habido buena convergencia: No se ha reducido el error entre el solido homogeneo y el metamaterial con esta ultima iteracion")
                else
                end
            end

        end
        #readline()

        hexaedros_analizados_metam_nuevo==[]
    end


    return hexaedros_analizados_metam_nuevo, inc_u_solido, inc_f_solido, u_solido_acumulado_nuevo, f_solido_acumulado_nuevo, fuerzas_hexaedros_acumulado_nuevo, energias_hexaedros_acumulado_nuevo, inc_energia_interna_total_solido, inc_energias_solidos_gruesos, energias_solidos_gruesos, inc_deform_equiv_solido_grueso, inc_u_metam, inc_f_metam, f_acumulada_elementos_nuevo, Youngs_tangentes_hexaedros_finos, Youngs_tangentes_hexaedros_gruesos, Youngs_tangentes_memoria_metam_nuevo, epsilon_voigt_nuevo, energias_acumuladas_elementos_nuevo, eps_x_ptos_int_nuevo, eps_xy_ptos_int_nuevo, eps_xz_ptos_int_nuevo, gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam, gdl_total_zona_hexaedros, gdl_fijos_zona_hexaedros, gdl_total_hexaedro, gdl_fijos_hexaedro, gdl_restring_totales_hexaedro_local, gdl_contorno_solido_completo_no_fronteras, energia_elemento_solido_grueso_metam, inc_energia_elemento_solido_grueso_metam, invariante_def_memoria_metam_nuevo
    
end #function resolucion_acoplamiento



function resolucion_iterativa_acoplamiento_zonas_mallado_fino_grueso_gradiente_fuerzas(hiperparametros ,X,Y,Z, N_3D_x_grueso, N_3D_y_grueso, N_3D_z_grueso, hexaedros_analizados_metam, hexaedros_mallado_fino_en_grueso, nodos_metam_en_hexaedros_finos, u_solido_acumulado, f_solido_acumulado, f_acumulado_hexaedros, energias_hexaedros_acumulado, epsilon_voigt_acumulado, f_acumulada_elementos_memoria, energias_acumuladas_elementos_memoria, eps_x_ptos_int_acumulado_memoria, eps_xy_ptos_int_acumulado_memoria, eps_xz_ptos_int_acumulado_memoria, coords_metam, conectividad_metam, coords_3D_fino, conectividad_3D_fino, conectividad_3D_grueso, volumenes_hexaedros_grueso, Youngs_elementos_finos_3D, Youngs_tangentes_hexaedros_gruesos, Youngs_tangentes_memoria_metam, tensiones_pandeo_elementos, punto_restriccion, gdl_restriccion, punto_desplazamiento, inc_valor_desplazamiento, coords_metam_hexaedros_exacto, conectiv_metam_hexaedros_exacto, porcentajes_barras_dentro_hexaedros, factor_barras_caras_solido, coords_metam_fronteras_exactas_hexaedros, gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam, gdl_no_fijos_total, gdl_total_zona_hexaedros, gdl_fijos_zona_hexaedros, gdl_no_fijos_zona_hexaedros, gdl_total_hexaedro, gdl_fijos_hexaedro, gdl_restring_fijos_hexaedro_local, gdl_restring_no_fijos_hexaedro_local, gdl_contorno_solido_completo_no_fronteras, invariante_def_memoria_metam, inc_valor_restricciones_no_fijas_total, K_step, gradiente_fuerzas_nodos, datos_contribucion_barras_gdl, iter_desplazamiento, radio)
    
    Young_solido_no_daño=curva_tension_deformacion_conjunto(10)

    # INICIALIZAR PARA EL ANALISIS DEL SOLIDO
    N_elementos_solidos_grueso=length(conectividad_3D_grueso[:,1])
    N_elementos_solidos_fino=length(conectividad_3D_fino[:,1])
    inc_u_solido=0; inc_f_solido=0; epsilon_voigt_nuevo=0; inc_energia_interna_total_solido=0
    young_material=Young_material(); SSR_antiguo=0
    hexaedros_analizados_metam_nuevo=copy(hexaedros_analizados_metam); # todos los elementos analizados en este desplazamiento
    inc_energia_hexaedros_total=0;
    energias_hexaedros_acumulado_nuevo=[]
    inc_deform_equiv_hexaed_aparente=zeros(N_elementos_solidos_grueso); inc_deform_equiv_solido_grueso=zeros(N_elementos_solidos_grueso); deform_equiv_solido_grueso=zeros(N_elementos_solidos_grueso) 
    inc_energias_solidos_gruesos=zeros(N_elementos_solidos_grueso); energias_solidos_gruesos=zeros(N_elementos_solidos_grueso)
    u_solido_acumulado_nuevo=[]; f_solido_acumulado_nuevo=[]; fuerzas_hexaedros_acumulado_nuevo=[]
    inc_energias_hexaedros_grueso_aparente=zeros(N_elementos_solidos_grueso)
    Youngs_tangentes_hexaedros_finos=copy(Youngs_elementos_finos_3D)
    
    
    # INICIALIZAR PARA EL ANALISIS DEL METAMATERIAL
    inc_u_metam=Dict(); inc_f_metam=Dict(); inc_energia_zona_metam=Dict()
    inc_energias_elementos=Dict(); inc_energia_elemento_solido_grueso_metam=Dict(); energia_elemento_solido_grueso_metam=Dict()
    Youngs_tangentes_memoria_metam_nuevo=Dict()
    f_acumulada_elementos_nuevo=Dict(); energias_acumuladas_elementos_nuevo=Dict()
    eps_x_ptos_int_nuevo=Dict(); eps_xy_ptos_int_nuevo=Dict(); eps_xz_ptos_int_nuevo=Dict()
    invariante_def_memoria_metam_nuevo=Dict()

    
    
    # HIPERPARAMETROS
    iter_max=hiperparametros[2];  tolerancia=hiperparametros[3]; tolerancia_2=hiperparametros[4]; diferencial=hiperparametros[5]
    eps_limite=hiperparametros[1]; modelo=hiperparametros[6]; iter_max_gradiente=hiperparametros[10]
    K_step_desplaz=hiperparametros[8]; K_step_giros=hiperparametros[9]

    ### si el exponente es cero, entonces, se les pone a todos los hexaedros finos el mismo Young que al hexaedro grueso
    exponente_reparto_young=hiperparametros[7]

    SSR_vector=[]; error_fuerzas=[]; gradiente_fuerzas_old=gradiente_fuerzas_nodos; factor=false; error_total=[]

    #K_step=copy(K_step_original)


    # COMO SE ANALIZAN TODOS LOS HEXAEDROS, DETERMINAMOS LOS GDL NO FIJOS TOTALES Y DE CADA HEXAEDRO
     
    
    ### GENERAR LOS GDL DE LAS ZONAS
    # SI NO SE HA ESTUDIADO TODAVIA EL HEXAEDRO (Y SU ZONA CORRESPONDIENTE ALREDEDOR, SE GENERAN LOS GDL CORRESPONDIENTES)

    # Si se cumple esta condicion, implica que es la primera vez que se calcula por gradiente
    if gdl_no_fijos_total==[]
        for elemento_solido_grueso in 1:N_elementos_solidos_grueso
            gdl_total_zona_hexaedros[elemento_solido_grueso], gdl_fijos_zona_hexaedros[elemento_solido_grueso], gdl_no_fijos_zona_hexaedros[elemento_solido_grueso] =clasificacion_gdl_zona_mallado_fino_grueso_gradiente(elemento_solido_grueso, coords_metam, conectividad_metam, conectiv_metam_hexaedros_exacto, coords_metam_hexaedros_exacto, coords_metam_fronteras_exactas_hexaedros, gdl_restring_globales_metam, gdl_contorno_solido_completo_no_fronteras)

            #println("TOTALES");println(gdl_total_zona_hexaedros[elemento_solido_grueso]);println()
            #println("fijos");println(gdl_fijos_zona_hexaedros[elemento_solido_grueso]);println()
            #println("no fijos");println(gdl_no_fijos_zona_hexaedros[elemento_solido_grueso]);println()
            #readline()

            # GDL CON RESTRICCIONES FIJAS DE LA ZONA 
            gdl_restring_fijos_hexaedro_local[elemento_solido_grueso]=indexin(gdl_fijos_zona_hexaedros[elemento_solido_grueso], gdl_total_zona_hexaedros[elemento_solido_grueso])

            # GDL CON RESTRICCIONES NO FIJAS DE LA ZONA 
            gdl_restring_no_fijos_hexaedro_local[elemento_solido_grueso]=indexin(gdl_no_fijos_zona_hexaedros[elemento_solido_grueso], gdl_total_zona_hexaedros[elemento_solido_grueso])


            # GDL NO FIJOS GLOBALES
            append!(gdl_no_fijos_total, gdl_no_fijos_zona_hexaedros[elemento_solido_grueso])
            unique!(sort!(gdl_no_fijos_total))
        end
    else
    end
    #println(gdl_no_fijos_total)


    # ITERAR PARA CONVERGER EN LA SOLUCION DEL MULTIESCALA
    
    for iter in 1:iter_max

        ######################### SOLIDO CONTINUO (MALLADO FINO) ######################################
        println("Iteracion desplazamiento ITER: ",iter ,"/", iter_max)
        println("inicio analisis del solido continuo (mallado fino)")

        ### TRADUCIR LAS CONDICIONES DE CONTORNO DEL SOLIDO A GRADOS DE LIBERTAD DEL SOLIDO
        gdl_restring_globales_solido, inc_valor_restricciones_solido=condic_contorno_solido(punto_restriccion, gdl_restriccion, punto_desplazamiento, inc_valor_desplazamiento, coords_3D_fino)

       
        ### RESOLVER LOS INCREMENTOS DE DESPLAZAMIENTOS DEL SOLIDO        
        inc_u_solido, inc_f_solido, inc_energias_hexaedros, deform_equiv_solido, inc_deform_equiv_solido, epsilon_voigt_nuevo = resolucion_solido_3D(coords_3D_fino, conectividad_3D_fino, Youngs_tangentes_hexaedros_finos, epsilon_voigt_acumulado, gdl_restring_globales_solido, inc_valor_restricciones_solido)
        
        ### INCREMENTO DE ENERGIA TOTAL DEL SOLIDO
        inc_energia_interna_total_solido=sum(inc_energias_hexaedros)        

        ### DESPLAZAMIENTOS Y FUERZAS TOTALES A PARTIR DE LOS INCREMENTOS 
        u_solido_acumulado_nuevo = u_solido_acumulado + inc_u_solido
        f_solido_acumulado_nuevo = f_solido_acumulado + inc_f_solido

        ### ENERGIA TOTAL ACUMULADA EN CADA ELEMENTO Y VECTOR DE FUERZAS EN CADA ELEMENTO
        energias_hexaedros_acumulado_nuevo, fuerzas_hexaedros_acumulado_nuevo = energias_y_fuerzas_acumuladas_elementos(energias_hexaedros_acumulado, inc_energias_hexaedros, f_acumulado_hexaedros, inc_u_solido, coords_3D_fino, conectividad_3D_fino, Youngs_tangentes_hexaedros_finos)
        
        #$$$
        #println();println("ITER desplazamiento: ",iter);
        println("deformacion equiv. max solido homogeneo= ",maximum(deform_equiv_solido));println("inc. energia solido homogeneo= ", sum(inc_energias_hexaedros))#;println("Youngs= ", Youngs_tangentes_hexaedros_finos);println()
        #println();println("Youngs_hexaed");println(Youngs_tangentes_hexaedros_finos);println(u_solido);println(fuerzas_nodos_solido)
        #readline()


        ############################ ENERGIAS Y DEFORMACIONES DE LOS SOLIDOS GRUESOS EN FUNCION DE LOS FINOS ##############################
        for elemento_solido_grueso in 1:N_elementos_solidos_grueso  

            inc_energias_solidos_gruesos[elemento_solido_grueso]=0
            energias_solidos_gruesos[elemento_solido_grueso]=0
            deform_equiv_solido_grueso[elemento_solido_grueso]=0; suma_2=0
            inc_deform_equiv_solido_grueso[elemento_solido_grueso]=0; suma_1=0

            N_hexaedros_finos_en_grueso=length(hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:])

            for elemento_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                inc_energias_solidos_gruesos[elemento_solido_grueso]=inc_energias_solidos_gruesos[elemento_solido_grueso] + inc_energias_hexaedros[elemento_fino]
                energias_solidos_gruesos[elemento_solido_grueso]=energias_solidos_gruesos[elemento_solido_grueso] + energias_hexaedros_acumulado_nuevo[elemento_fino]
                suma_1=suma_1 + inc_deform_equiv_solido[elemento_fino]^2/N_hexaedros_finos_en_grueso
                suma_2=suma_2 + deform_equiv_solido[elemento_fino]^2/N_hexaedros_finos_en_grueso
            end
            inc_deform_equiv_solid=sqrt(suma_1)#; println("1 ",sqrt(suma_1))
            deform_equiv_solido_grueso[elemento_solido_grueso]=sqrt(suma_2)
            inc_deform_equiv_solido_grueso[elemento_solido_grueso]=sqrt(inc_energias_solidos_gruesos[elemento_solido_grueso]/(1/2*Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]*volumenes_hexaedros_grueso[elemento_solido_grueso]))
            #println("2 ", inc_deform_equiv_solido_grueso[elemento_solido_grueso])
        end
        #readline()

        ################################ MULTIESCALA ####################################

        ### PASAR LAS RESTRICCIONES DEL SOLIDO A RESTRICCIONES DEL METAMATERIAL (CON 6 GDL EN VEZ DE 3)

        # if gdl_restring_globales_metam==[]
        #     N_restricciones = length(gdl_restriccion[:,1])
        #     N_desplazamientos = length(inc_valor_desplazamiento[:,1])

        #     punto_restriccion_metam = punto_restriccion
        #     gdl_restriccion_metam = [gdl_restriccion zeros(N_restricciones,3)]
        #     punto_desplazamiento_metam = punto_desplazamiento
        #     inc_valor_desplazamiento_metam = [inc_valor_desplazamiento zeros(N_desplazamientos,3)]

        #     # GRADOS DE LIBERTAD CON ALGUNA RESTRICCION Y EL VALOR DE ESA RESTRICCION (PUEDEN SER DEL CONTORNO O NO) del solido completo
        #     gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam=condic_contorno(punto_restriccion_metam, gdl_restriccion_metam ,punto_desplazamiento_metam ,inc_valor_desplazamiento_metam ,coords_metam)
        #     #println(gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam)

            
        #     # OBTENER LOS GDL DEL CONTORNO COMPLETO QUE NO ESTEN EN LAS FRONTERAS EN CONTACTO CON OTRAS ZONAS
        #     Nodos_metam_contorno_solido, gdl_contorno_solido_completo_no_fronteras=gdl_metamaterial_contorno_solido_no_fronteras(X,Y,Z,N_3D_x_grueso, N_3D_y_grueso, N_3D_z_grueso, coords_metam)

        # else
        # end


        ### OBTENER LAS CONDICIONES DE CONTORNO PARA LOS NODOS RESTRINGIDOS NO FIJOS UTILIZANDO TANTO LOS DESPLAZAMIENTOS DEL SOLIDO COMO LAS RESTRICCIONES GLOBALES
        if inc_valor_restricciones_no_fijas_total==[]
            N_gdl_restricciones_no_fijas_totales=length(gdl_no_fijos_total)
            inc_valor_restricciones_no_fijas_total=zeros(N_gdl_restricciones_no_fijas_totales)
            K_step=zeros(N_gdl_restricciones_no_fijas_totales)
            factor=true

            inc_u_solido_condic_contorno_multiescala=inc_u_solido

            for i in 1:N_gdl_restricciones_no_fijas_totales
                gdl=gdl_no_fijos_total[i]
                #println(gdl)

                # DETERMINAR A QUE NODO LE CORRESPONDE CADA GDL Y QUE VALOR DE DESPLAZAMIENTO LE TOCA
                nodo=Int64(ceil(gdl/6))

                for elemento_solido_fino in 1:N_elementos_solidos_fino
                    indice=findall(in(nodo), nodos_metam_en_hexaedros_finos[elemento_solido_fino])
                    
                    # SI ESTÁ EN EL ELEMENTO SOLIDO QUE ESTAMOS MIRANDO AHORA
                    if indice!=[]
                        indice_nodo=findall(in(nodo),coords_metam[:,4])
                        coords_punto=coords_metam[indice_nodo[1],1:3]'
                        
                        resto=gdl%6 # como imponemos desplazamientos y giros, el resto irá entra 0 y 5

                        # SI EL GDL CORRESPONDE A UN DESPLAZAMIENTO EN X
                        if resto == 1
                            inc_desplaz_nodo=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                            k_step=K_step_desplaz

                        # SI EL GDL CORRESPONDE A UN DESPLAZAMIENTO EN Y    
                        elseif resto == 2
                            inc_desplaz_nodo=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                            k_step=K_step_desplaz

                        # SI EL GDL CORRESPONDE A UN DESPLAZAMIENTO EN Z    
                        elseif resto == 3
                            inc_desplaz_nodo=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                            k_step=K_step_desplaz

                        # SI EL GDL CORRESPONDE A UN GIRO EN X (AROXIMAMOS EL GIRO)   
                        elseif resto == 4                        
                            inc_desplaz_nodo_y=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                            inc_desplaz_nodo_diferencial_y=desplazamientos_elemento(coords_punto+[0 0 diferencial], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                            inc_desplaz_nodo_1=-(inc_desplaz_nodo_diferencial_y-inc_desplaz_nodo_y)/ diferencial

                            inc_desplaz_nodo_z=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                            inc_desplaz_nodo_diferencial_z=desplazamientos_elemento(coords_punto+[0 diferencial 0], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                            inc_desplaz_nodo_2=(inc_desplaz_nodo_diferencial_z-inc_desplaz_nodo_z)/ diferencial
        
                            inc_desplaz_nodo=(inc_desplaz_nodo_1+inc_desplaz_nodo_2)/2
                            k_step=K_step_giros

                        # SI EL GDL CORRESPONDE A UN GIRO EN Y (AROXIMAMOS EL GIRO)   
                        elseif resto == 5                        
                            inc_desplaz_nodo_z=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                            inc_desplaz_nodo_diferencial_z=desplazamientos_elemento(coords_punto+[diferencial 0 0], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                            inc_desplaz_nodo_1=-(inc_desplaz_nodo_diferencial_z-inc_desplaz_nodo_z)/ diferencial
                            
                            inc_desplaz_nodo_x=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                            inc_desplaz_nodo_diferencial_x=desplazamientos_elemento(coords_punto+[0 0 diferencial], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                            inc_desplaz_nodo_2=(inc_desplaz_nodo_diferencial_x-inc_desplaz_nodo_x)/ diferencial
                            
                            inc_desplaz_nodo=(inc_desplaz_nodo_1+inc_desplaz_nodo_2)/2
                            k_step=K_step_giros

                        # SI EL GDL CORRESPONDE A UN GIRO EN Z (AROXIMAMOS EL GIRO)   
                        elseif resto == 0                     
                            inc_desplaz_nodo_y=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                            inc_desplaz_nodo_diferencial_y=desplazamientos_elemento(coords_punto+[diferencial 0 0], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                            inc_desplaz_nodo_1=(inc_desplaz_nodo_diferencial_y-inc_desplaz_nodo_y)/ diferencial
                            
                            inc_desplaz_nodo_x=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                            inc_desplaz_nodo_diferencial_x=desplazamientos_elemento(coords_punto+[0 diferencial 0 ], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                            inc_desplaz_nodo_2=-(inc_desplaz_nodo_diferencial_x-inc_desplaz_nodo_x)/ diferencial

                            inc_desplaz_nodo=(inc_desplaz_nodo_1+inc_desplaz_nodo_2)/2
                            k_step=K_step_giros
                        end

                        inc_valor_restricciones_no_fijas_total[i]=inc_desplaz_nodo
                        K_step[i]=k_step
                        #K_step_original=K_step
                        #println(gdl," ",inc_desplaz_nodo," ",elemento_solido_fino) 
                        break
                    else
                    end
                end
            end
        else
        end
        #readline()

        for iter_gradiente in 1:iter_max_gradiente

            println();println("iter gradiente: ", Int64(iter_gradiente),"/",Int64(iter_max_gradiente)); println()

            ### SI EN ALGUN ELEMENTO SOLIDO DEL MALLADO FINO SE SOBREPASA LA TENSION LIMITE (ARBITRARIA) SE HACE EL ANALISIS MULTIESCALA DEL ELEMENTO DEL MALLADO GRUESO QUE CONTENGA A ESE ELEMENTO DEL MALLADO FINO
            for elemento_solido_grueso in 1:N_elementos_solidos_grueso  
        
                println();println("ELEMENTO ",elemento_solido_grueso, " Y SU ZONA")#;println("deform_equiv_tet= ", deform_equiv_solido[elemento_solido_grueso]);println()
                    

                # SI TODAVIA NO SE HA ESTUDIADO ESTA ZONA, SE PONE LA CONDICION DE CONTORNO DEL DESPLAZAMIENTO TOTAL DEL SOLIDO
                if findall(in(elemento_solido_grueso), hexaedros_analizados_metam)==[] #si no ha sido estudiada esta zona para desplazamientos iteraciones anteriores
                    println("nuevo elemento estudiado con multiescala en este desplazamiento")
                    inc_u_solido_condic_contorno_multiescala=u_solido_acumulado_nuevo
                    inc_valor_restricciones_fijos_global_metam_aparente=copy(inc_valor_restricciones_fijos_global_metam)*iter_desplazamiento

                # SI SE HA ESTUDIADO YA ESTA ZONA, LE IMPONEMOS EL INCREMENTO DE DESPLAZAMIENTO DEL SOLIDO
                else
                    inc_u_solido_condic_contorno_multiescala=inc_u_solido
                    inc_valor_restricciones_fijos_global_metam_aparente=copy(inc_valor_restricciones_fijos_global_metam)
                end


                ### OBTENER LAS CONDICIONES DE CONTORNO PARA LOS NODOS RESTRINGIDOS FIJOS UTILIZANDO TANTO LOS DESPLAZAMIENTOS DEL SOLIDO COMO LAS RESTRICCIONES GLOBALES
                N_gdl_restricciones_fijas_zona=length(gdl_fijos_zona_hexaedros[elemento_solido_grueso])
                inc_valor_restricciones_fijas_zona=zeros(N_gdl_restricciones_fijas_zona)

                for i in 1:N_gdl_restricciones_fijas_zona
                    gdl=gdl_fijos_zona_hexaedros[elemento_solido_grueso][i]
                    #println(gdl)

                    # Si el gdl esta en las restricciones globales
                    if (indice=findall(in(gdl),gdl_restring_globales_metam))!=[]
                        inc_valor_restricciones_fijas_zona[i]=inc_valor_restricciones_fijos_global_metam_aparente[indice[1]]

                    # Si el gdl no esta en las restricciones globales, se le restringe a partir de los desplazamientos del solido
                    else
                        # DETERMINAR A QUE NODO LE CORRESPONDE CADA GDL Y QUE VALOR DE DESPLAZAMIENTO LE TOCA
                        nodo=Int64(ceil(gdl/6))

                        for elemento_solido_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                            indice=findall(in(nodo), nodos_metam_en_hexaedros_finos[elemento_solido_fino])
                            
                            # SI ESTÁ EN EL ELEMENTO SOLIDO DE LA ZONA QUE ESTAMOS MIRANDO AHORA
                            if indice!=[]
                                indice_nodo=findall(in(nodo),coords_metam_fronteras_exactas_hexaedros[elemento_solido_grueso][:,4])
                                coords_punto=coords_metam_fronteras_exactas_hexaedros[elemento_solido_grueso][indice_nodo[1],1:3]'
                                
                                resto=gdl%6 # como imponemos desplazamientos y giros, el resto irá entra 0 y 5

                                # SI EL GDL CORRESPONDE A UN DESPLAZAMIENTO EN X
                                if resto == 1
                                    inc_desplaz_nodo=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                                    
                                # SI EL GDL CORRESPONDE A UN DESPLAZAMIENTO EN Y    
                                elseif resto == 2
                                    inc_desplaz_nodo=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                                    
                                # SI EL GDL CORRESPONDE A UN DESPLAZAMIENTO EN Z    
                                elseif resto == 3
                                    inc_desplaz_nodo=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                                    
                                # SI EL GDL CORRESPONDE A UN GIRO EN X (AROXIMAMOS EL GIRO)   
                                elseif resto == 4                        
                                    inc_desplaz_nodo_y=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                                    inc_desplaz_nodo_diferencial_y=desplazamientos_elemento(coords_punto+[0 0 diferencial], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                                    inc_desplaz_nodo_1=-(inc_desplaz_nodo_diferencial_y-inc_desplaz_nodo_y)/ diferencial

                                    inc_desplaz_nodo_z=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                                    inc_desplaz_nodo_diferencial_z=desplazamientos_elemento(coords_punto+[0 diferencial 0], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                                    inc_desplaz_nodo_2=(inc_desplaz_nodo_diferencial_z-inc_desplaz_nodo_z)/ diferencial
            
                                    inc_desplaz_nodo=(inc_desplaz_nodo_1+inc_desplaz_nodo_2)/2
                                    
                                # SI EL GDL CORRESPONDE A UN GIRO EN Y (AROXIMAMOS EL GIRO)   
                                elseif resto == 5                        
                                    inc_desplaz_nodo_z=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                                    inc_desplaz_nodo_diferencial_z=desplazamientos_elemento(coords_punto+[diferencial 0 0], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[3]
                                    inc_desplaz_nodo_1=-(inc_desplaz_nodo_diferencial_z-inc_desplaz_nodo_z)/ diferencial
                                    
                                    inc_desplaz_nodo_x=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                                    inc_desplaz_nodo_diferencial_x=desplazamientos_elemento(coords_punto+[0 0 diferencial], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                                    inc_desplaz_nodo_2=(inc_desplaz_nodo_diferencial_x-inc_desplaz_nodo_x)/ diferencial
                                    
                                    inc_desplaz_nodo=(inc_desplaz_nodo_1+inc_desplaz_nodo_2)/2

                                # SI EL GDL CORRESPONDE A UN GIRO EN Z (AROXIMAMOS EL GIRO)   
                                elseif resto == 0                     
                                    inc_desplaz_nodo_y=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                                    inc_desplaz_nodo_diferencial_y=desplazamientos_elemento(coords_punto+[diferencial 0 0], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[2]
                                    inc_desplaz_nodo_1=(inc_desplaz_nodo_diferencial_y-inc_desplaz_nodo_y)/ diferencial
                                    
                                    inc_desplaz_nodo_x=desplazamientos_elemento(coords_punto, elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                                    inc_desplaz_nodo_diferencial_x=desplazamientos_elemento(coords_punto+[0 diferencial 0 ], elemento_solido_fino, coords_3D_fino, conectividad_3D_fino, inc_u_solido_condic_contorno_multiescala)[1]
                                    inc_desplaz_nodo_2=-(inc_desplaz_nodo_diferencial_x-inc_desplaz_nodo_x)/ diferencial

                                    inc_desplaz_nodo=(inc_desplaz_nodo_1+inc_desplaz_nodo_2)/2
                                end

                                inc_valor_restricciones_fijas_zona[i]=inc_desplaz_nodo
                                #println(gdl," ",inc_desplaz_nodo," ",elemento_solido_fino) 
                                
                                break
                            else
                            end
                        end
                    end
                end

                
                ### OBTENER LAS CONDICIONES DE CONTORNO PARA LOS NODOS RESTRINGIDOS NO FIJOS UTILIZANDO TANTO LOS DESPLAZAMIENTOS DEL SOLIDO COMO LAS RESTRICCIONES GLOBALES
                N_gdl_restricciones_no_fijas_zona=length(gdl_no_fijos_zona_hexaedros[elemento_solido_grueso])
                inc_valor_restricciones_no_fijas_zona=zeros(N_gdl_restricciones_no_fijas_zona)

                for i in 1:N_gdl_restricciones_no_fijas_zona
                    gdl=gdl_no_fijos_zona_hexaedros[elemento_solido_grueso][i]
                    #println(gdl)

                    # DETERMINAR A QUE NODO LE CORRESPONDE CADA GDL Y QUE VALOR DE DESPLAZAMIENTO LE TOCA
                    indice=findall(in(gdl), gdl_no_fijos_total);# println(indice)

                    inc_valor_restricciones_no_fijas_zona[i]=inc_valor_restricciones_no_fijas_total[indice[1]]
                    #println(gdl," ",inc_desplaz_nodo," ",elemento_solido_fino) 

                end


                ### OBTENER TODAS LAS CONDICIONES DE CONTORNO DE TODOS LOS GDL RESTRINGIDOS
                gdl_restring_totales_hexaedro_local=vcat(gdl_restring_fijos_hexaedro_local[elemento_solido_grueso], gdl_restring_no_fijos_hexaedro_local[elemento_solido_grueso])
                inc_valor_restricciones_totales_zona=vcat(inc_valor_restricciones_fijas_zona, inc_valor_restricciones_no_fijas_zona)



                ################################# RESOLUCION DE LA ZONA ##################################

                ########### OBTENER LOS GDL QUE VAN A INTERVENIR Y LA CONECTIVIDAD DENTRO DEL hexaedro
                # coordenadas de los nodos de la zona con metamaterial (un elemento hexaedrico)
                coords_metam_elemento=coords_metam_hexaedros_exacto[elemento_solido_grueso]

                # Conectividad de las barras de la zona
                conectividad_zona_metam=conectiv_metam_hexaedros_exacto[elemento_solido_grueso]
                #porcentajes_barras_zona=porcentajes_barras_dentro_hexaedros[elemento_solido_zona]
                N_elementos_metam=length(conectividad_zona_metam[:,1])
                #$$$$
                #println();println("nodos_tot_extend");println(coords_metam_elemento_extendido);println()


                ########### INICIALIZAR VARIABLES Y DECIDIR QUE CONDICIONES DE CONTORNO
                # INICIALIZAR LA MATRIZ DE RIGIDEZ DE LA ZONA DE METAMATERIAL QUE ESTAMOS ESTUDIANDO

                # SI NO HA SIDO ESTUDIADA ESTA ZONA PARA DESPLAZAMIENTOS ITERACIONES ANTERIORES
                if findall(in(elemento_solido_grueso), hexaedros_analizados_metam)==[] 
                    
                    # si no hay registros anteriores, inicializamos las matrices
                    inicializar_cero=Dict(); eps_x_ptos_int_acumulado=Dict(); eps_xy_ptos_int_acumulado=Dict(); eps_xz_ptos_int_acumulado=Dict(); f_acumulada_elementos_zona=Dict() ;energias_acumuladas_elementos_zona=Dict(); invariante_def_metam=Dict()
                    for i in 1:N_elementos_metam
                        element_metam=conectividad_zona_metam[i,3]
                        inicializar_cero[element_metam]=0 # un diccionario inicializado a cero para uso generico
                        eps_x_ptos_int_acumulado[element_metam]=[0,0]      
                        eps_xy_ptos_int_acumulado[element_metam]=[0,0]
                        eps_xz_ptos_int_acumulado[element_metam]=[0,0]
                        energias_acumuladas_elementos_zona[element_metam]=0 
                        f_acumulada_elementos_zona[element_metam] =zeros(12)
                        invariante_def_metam[element_metam]=0    
                    end
                    Youngs_tangentes_in=curva_tension_deformacion_tabla_Osgood(inicializar_cero,"deform_equiv",inicializar_cero)  #si no se tienen registros de modulos de young anteriores, se empieza por el valor base del material
                    Youngs_tangentes_memoria_metam[elemento_solido_grueso]=Youngs_tangentes_in
                    inc_deform_equiv_hexaed_aparente[elemento_solido_grueso]=deform_equiv_solido_grueso[elemento_solido_grueso]
                    inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]=energias_solidos_gruesos[elemento_solido_grueso]
                    
                # SI LA ZONA YA HA SIDO ESTUDIADA, COGEMOS LOS VALORES QUE YA EXISTEN Y ESTUDIAMOS INCREMENTOS DESDE ESE PUNTO
                else                       
                    #si ya hay un vector de Youngs calculado para el desplazamiento anterior, se coje ese para facilitar la convergencia (cambio menos brusco)          
                    Youngs_tangentes_in=Dict(); eps_x_ptos_int_acumulado=Dict(); eps_xy_ptos_int_acumulado=Dict(); eps_xz_ptos_int_acumulado=Dict(); f_acumulada_elementos_zona=Dict() ;energias_acumuladas_elementos_zona=Dict(); invariante_def_metam=Dict()
                    
                    Youngs_tangentes_in=copy(Youngs_tangentes_memoria_metam[elemento_solido_grueso])
                    eps_x_ptos_int_acumulado=copy(eps_x_ptos_int_acumulado_memoria[elemento_solido_grueso])      
                    eps_xy_ptos_int_acumulado=copy(eps_xy_ptos_int_acumulado_memoria[elemento_solido_grueso])   
                    eps_xz_ptos_int_acumulado=copy(eps_xz_ptos_int_acumulado_memoria[elemento_solido_grueso])
                    energias_acumuladas_elementos_zona =copy(energias_acumuladas_elementos_memoria[elemento_solido_grueso])
                    f_acumulada_elementos_zona =copy(f_acumulada_elementos_memoria[elemento_solido_grueso])
                    invariante_def_metam=copy(invariante_def_memoria_metam[elemento_solido_grueso])

                    inc_deform_equiv_hexaed_aparente[elemento_solido_grueso]=inc_deform_equiv_solido_grueso[elemento_solido_grueso]
                    inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]=inc_energias_solidos_gruesos[elemento_solido_grueso]
                end
                                    



                ########### RESOLVER EL PROBLEMA DEL HEXAEDRO CON METAMATERIAL (EL ZOOM EN ESE ELEMENTO DEL SOLIDO)                 

                # obtener los incrementos en desplazamientos, la energia de la zona extedida y los modulos de Young
                inc_u_metam[elemento_solido_grueso], inc_f_metam[elemento_solido_grueso], f_acumulada_elementos_nuevo[elemento_solido_grueso], inc_energia_zona_metam[elemento_solido_grueso], inc_energias_elementos[elemento_solido_grueso], energias_acumuladas_elementos_nuevo[elemento_solido_grueso], Youngs_tangentes_memoria_metam_nuevo[elemento_solido_grueso], eps_x_ptos_int_nuevo[elemento_solido_grueso], eps_xy_ptos_int_nuevo[elemento_solido_grueso], eps_xz_ptos_int_nuevo[elemento_solido_grueso], invariante_def_memoria_metam_nuevo[elemento_solido_grueso]=resolucion_iterativa_metamaterial(coords_metam_elemento, conectividad_zona_metam, gdl_restring_totales_hexaedro_local, inc_valor_restricciones_totales_zona , Youngs_tangentes_in, tensiones_pandeo_elementos[elemento_solido_grueso], eps_x_ptos_int_acumulado, eps_xy_ptos_int_acumulado, eps_xz_ptos_int_acumulado, energias_acumuladas_elementos_zona, f_acumulada_elementos_zona, invariante_def_metam, radio)
                # println("desplazamientos y fuerzas")
                # for i in 1:Int64(length(inc_f_metam[elemento_solido_grueso])/6)
                #     println(i," ",inc_u_metam[elemento_solido_grueso][6*i-5:6*i]," ",inc_f_metam[elemento_solido_grueso, elemento_solido_grueso_zona][6*i-5:6*i])
                # end
                # println()
                # println("energia")
                # for i in conectividad_zona_metam[:,3]
                #     #println(i," ",inc_energias_elementos[elemento_solido_grueso][i], " ", energias_acumuladas_elementos_nuevo[elemento_solido_grueso][i], " ", Youngs_tangentes_memoria_metam_nuevo[elemento_solido_grueso][i])
                #     indice=findall(in(i),conectividad_zona_metam[:,3])
                #     nodos=conectividad_zona_metam[indice[1],1:2]
                #     println(nodos)
                #     indices=findall(in(nodos),coords_metam_hexaedros_exacto[elemento_solido_grueso][:,4])
                #     println(f_acumulada_elementos_nuevo[elemento_solido_grueso][i]-f_acumulada_elementos_memoria[elemento_solido_grueso][i])
                #     println("   ",inc_u_metam[elemento_solido_grueso][6*indices[1]-5:6*indices[1]])
                #     println("   ",inc_u_metam[elemento_solido_grueso][6*indices[2]-5:6*indices[2]])
                # end

                # readline()


                ######## UNA VEZ ANALIZADOS TODOS LOS ELEMENTOS DE LA ZONA, SOLO NOS IMPORTA LA ENERGIA DEL ELEMENTO SOLIDO GRUESO

                # sumamos la energia de las barras del elemento
                inc_energia_hexaedro_metam=0; energia_hexaedro_metam=0; inc_energia_daño_cero=0
                for barra in conectiv_metam_hexaedros_exacto[elemento_solido_grueso][:,3]
                    inc_energia_hexaedro_metam=inc_energia_hexaedro_metam + inc_energias_elementos[elemento_solido_grueso][barra]*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]
                    energia_hexaedro_metam=energia_hexaedro_metam + energias_acumuladas_elementos_nuevo[elemento_solido_grueso][barra]*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]
                    
                    inc_energia_daño_cero=inc_energia_daño_cero+inc_energias_elementos[elemento_solido_grueso][barra]*young_material/Youngs_tangentes_memoria_metam_nuevo[elemento_solido_grueso][barra]*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]
                    #$$$$$
                    # indice=findall(in(barra),conectividad_metam[:,3])
                    # nodo_1=conectividad_metam[indice,1]; nodo_2=conectividad_metam[indice,2]
                    # println(barra," ", nodo_1," ", nodo_2, " ", porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]," ",inc_energias_elementos[elemento_solido_grueso][barra]*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra])
                end
                daño_hexaedro_grueso=inc_energia_hexaedro_metam/inc_energia_daño_cero

                inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]=inc_energia_hexaedro_metam
                energia_elemento_solido_grueso_metam[elemento_solido_grueso]=energia_hexaedro_metam
                
                println();println("hexaedro= ",elemento_solido_grueso," inc_energia_exacto= ",inc_energia_hexaedro_metam)
                #readline()


                ########### OBTENER EL PARAMETRO DE DAÑO QUE HEMOS DEFINIDO PARA EL REPARTO DEL YOUNG A CADA ELEMENTO FINO

                ######## PARAMETRO DE DAÑO
                if (exponente_reparto_young==0) && (modelo != "modelo_daño_hexaedros_finos")
                else
                    daño_hexaedros_finos=Dict()
                    
                    for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                        inc_energia_hexaedros_finos=0; inc_energia_daño_nulo=0; barras=[]; daño_hexaedros_finos[hexaedro_fino]=0

                        for nodo in nodos_metam_en_hexaedros_finos[hexaedro_fino]
                            indices=findall(in(nodo), conectividad_zona_metam[:,1:2])
                            for indice in indices
                                push!(barras, conectividad_zona_metam[indice[1],3])
                            end
                            unique!(sort!(barras))
                        end
                        
                        for barra in barras                    
                            inc_energia_hexaedros_finos=inc_energia_hexaedros_finos + inc_energias_elementos[elemento_solido_grueso][barra]#*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]
                            inc_energia_daño_nulo=inc_energia_daño_nulo + inc_energias_elementos[elemento_solido_grueso][barra]*young_material/Youngs_tangentes_memoria_metam_nuevo[elemento_solido_grueso][barra]#*porcentajes_barras_dentro_hexaedros[elemento_solido_grueso][barra]
                        end

                        daño_hexaedros_finos[hexaedro_fino]=inc_energia_hexaedros_finos/inc_energia_daño_nulo
                        #println(daño_hexaedros_finos[hexaedro_fino])
                    end
                end


                ########### OBTENER EL MODULO DE YOUNG TANGENTE NUEVO DE ESTE HEXAEDRO (APROXIMACION) --> IGUALAMOS LAS ENERGIAS, Y SUPONIENDO LOS DESPLAZAMIENTOS BIEN CALCULADOS, OBTENEMOS EL NUEVO YOUNG
                # Lo que hacemos es como calcular el modulo de Young tangente del hexaedro pero a partir de igualar la energia

                ### DEPENDIENDO DEL MODELO UTILIZADO SE ACTUALIZAN LOS YOUNGS TANGENTES DE UNA FORMA U OTRA
                if modelo == "modelo_energia_hexaedros_gruesos"
                    # Young_tangente_nuevo_2=inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]/(1/2*inc_deform_equiv_hexaed_aparente[elemento_solido_grueso]^2*volumenes_hexaedros_grueso[elemento_solido_grueso]) # nueva estimacion del valor del Young del hexaedro
                    Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]=inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]/(inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]/Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso])
                    
                    ### si el exponente es cero, entonces, se les pone a todos los hexaedros finos el mismo Young que al hexaedro grueso
                    if exponente_reparto_young==0
                        for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                            Youngs_tangentes_hexaedros_finos[hexaedro_fino]=Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]         
                        end
                    else
                        for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                            Youngs_tangentes_hexaedros_finos[hexaedro_fino]=Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]*(daño_hexaedros_finos[hexaedro_fino]/daño_hexaedro_grueso)^exponente_reparto_young          
                        end
                    end

                    println("Young hexaedro grueso modelo energia: ",Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso], " ", sum(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])/length(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]]))
                    println("Young hexaedro grueso modelo daño: ",Young_solido_no_daño*daño_hexaedro_grueso, " ", sum(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])/length(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]]))
                    println("youngs hexaedros finos: ", Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])
                    ##Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]].= Young_tangente_nuevo# (Youngs_elementos_finos_3D[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,1]] + Young_tangente_nuevo)/2
                
                elseif modelo == "modelo_daño_hexaedros_finos"
                 
                    for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                        Youngs_tangentes_hexaedros_finos[hexaedro_fino]=Young_solido_no_daño*daño_hexaedros_finos[hexaedro_fino]             
                    end
                                       
                    Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]=Young_solido_no_daño*daño_hexaedro_grueso
                    println("Young hexaedro grueso modelo energia: ",inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]/(inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]/Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]))
                    println("Young hexaedro grueso modelo daño: ",Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso], " ", sum(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])/length(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]]))
                    println("youngs hexaedros finos: ", Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])
                
                elseif modelo == "modelo_daño_hexaedros_gruesos"

                    if exponente_reparto_young==0
                        for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                            Youngs_tangentes_hexaedros_finos[hexaedro_fino]=Young_solido_no_daño*daño_hexaedro_grueso
                        end
                    else
                        for hexaedro_fino in hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]
                            Youngs_tangentes_hexaedros_finos[hexaedro_fino]=Young_solido_no_daño*daño_hexaedro_grueso*(daño_hexaedros_finos[hexaedro_fino]/daño_hexaedro_grueso)^exponente_reparto_young   
                        end
                    end
                    
                    Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]=Young_solido_no_daño*daño_hexaedro_grueso
                    println("Young hexaedro grueso modelo energia: ",inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]/(inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]/Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso]))
                    println("Young hexaedro grueso modelo daño: ",Youngs_tangentes_hexaedros_gruesos[elemento_solido_grueso], " ", sum(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])/length(Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]]))
                    println("youngs hexaedros finos: ", Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,:]])
                
                else
                    println("No se ha introducido modelo de daño")
                end
                    #readline()
                #println(1/2*inc_deform_equiv_hexaed_aparente[elemento_solido_grueso]^2*volumenes_hexaedros_grueso[elemento_solido_grueso]," ", inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]/Youngs_tangentes_hexaedros_finos[hexaedros_mallado_fino_en_grueso[elemento_solido_grueso,1]])
                #println(Young_tangente_nuevo," ",Young_tangente_nuevo_2," ",Youngs_elementos_finos_3D[elemento_solido_grueso])    
                
                #$$$$
                println("inc_energ_tet= ",inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso]);println("inc_energ_met= ",inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]);println()
                println("energ_tet= ",energias_solidos_gruesos[elemento_solido_grueso]);println("energ_met= ",energia_elemento_solido_grueso_metam[elemento_solido_grueso]);println()
                #readline()           
                
            end
        
            println(Youngs_tangentes_hexaedros_gruesos)
            unique!(sort!(hexaedros_analizados_metam_nuevo))


            ### METODO GRADIENTE: ACTUALIZAR INC DESPLAZAMIENTOS NO FIJOS: Si en este fichero se introduce en "HIPERPARAMETROS": gradiente=true, en cada step de desplazamientos, aplica el metodo del gradiente para corregir que las fuerzas deben ser continuas en dos hexaedros adyacentes (no solo los desplazamientos) --> la funcion objetivo a optimizar es la que se consigue sumando las fuerzas en los nodos comunes que estan en dos hexaedros adyacentes
            gradiente_fuerzas_nodos=zeros(length(gdl_no_fijos_total)) #; fuerzas_gdl=zeros(length(gdl_no_fijos_total))


            if collect(keys(datos_contribucion_barras_gdl))==[]
                if length(hexaedros_analizados_metam_nuevo)==N_elementos_solidos_grueso
                    for hexaedro_grueso in 1:N_elementos_solidos_grueso
                        for gdl in gdl_no_fijos_zona_hexaedros[hexaedro_grueso]

                            # if (gdl%6==1) || (gdl%6==2) || (gdl%6==3) 
                            #     K_step=1
                            # else
                            #     K_step=1#/(10^2)
                            # end

                            nodo=ceil(gdl/6)
                            indice=findall(in(gdl),gdl_no_fijos_total)[1]
                            
                            barras_nodo_local=findall(in(nodo), conectiv_metam_hexaedros_exacto[hexaedro_grueso][:,1:2])


                            for barra_local in barras_nodo_local 
                                gdl_local_barra=gdl%6 
                                if gdl_local_barra!=0
                                    gdl_local_barra=gdl_local_barra+6*(barra_local[2]-1)
                                else
                                    gdl_local_barra=6+6*(barra_local[2]-1)
                                end
                                barra=conectiv_metam_hexaedros_exacto[hexaedro_grueso][barra_local[1],3]
 
                                gradiente_fuerzas_nodos[indice]=gradiente_fuerzas_nodos[indice]+(f_acumulada_elementos_nuevo[hexaedro_grueso][barra][gdl_local_barra]-f_acumulada_elementos_memoria[hexaedro_grueso][barra][gdl_local_barra])*porcentajes_barras_dentro_hexaedros[hexaedro_grueso][barra]*factor_barras_caras_solido[barra]
                                
                                fila_con_datos=[indice hexaedro_grueso barra gdl_local_barra]
                                if (haskey(datos_contribucion_barras_gdl, gdl)==false)
                                    datos_contribucion_barras_gdl[gdl]=fila_con_datos
                                else
                                    datos_contribucion_barras_gdl[gdl]=vcat(datos_contribucion_barras_gdl[gdl], fila_con_datos)
                                end
 
                            end
                        end
                    end
                else
                end
            else
                for gdl in gdl_no_fijos_total
                    for i in 1:length(datos_contribucion_barras_gdl[gdl][:,1])
                        indice=datos_contribucion_barras_gdl[gdl][i,1]
                        hexaedro_grueso=datos_contribucion_barras_gdl[gdl][i,2]
                        barra=datos_contribucion_barras_gdl[gdl][i,3]
                        gdl_local_barra=datos_contribucion_barras_gdl[gdl][i,4]
                        gradiente_fuerzas_nodos[indice]=gradiente_fuerzas_nodos[indice]+(f_acumulada_elementos_nuevo[hexaedro_grueso][barra][gdl_local_barra]-f_acumulada_elementos_memoria[hexaedro_grueso][barra][gdl_local_barra])*porcentajes_barras_dentro_hexaedros[hexaedro_grueso][barra]*factor_barras_caras_solido[barra]
                        
                    end
                end
            end
            
            ####
            #inc_valor_restricciones_no_fijas_total_old=inc_valor_restricciones_no_fijas_total
            #inc_valor_restricciones_no_fijas_total=inc_valor_restricciones_no_fijas_total-0.015*gradiente_fuerzas_nodos.*(abs.(inc_valor_restricciones_no_fijas_total)./fuerzas_gdl)
            
            #inc_valor_restricciones_no_fijas_total=inc_valor_restricciones_no_fijas_total-step_gradiente*gradiente_fuerzas_nodos
            
            
            #plotlyjs()
            mejora_error_inferior=0;mejora_error_superior=0
            if factor==true      
                #println(K_step)   
                inc_valor_restricciones_no_fijas_total=inc_valor_restricciones_no_fijas_total-K_step.*gradiente_fuerzas_nodos 
                factor=false
            else     
                reduccion_error_fuerzas= gradiente_fuerzas_nodos./gradiente_fuerzas_old*100
                mejora_error_inferior=(abs.(reduccion_error_fuerzas).>60)
                mejora_error_superior=(abs.(reduccion_error_fuerzas).<200)#300
                mejora_error_negativa=(reduccion_error_fuerzas.<0)#-30
                #println(mejora_error_inferior)
                K_step=(.!mejora_error_inferior .*K_step+ 3*mejora_error_inferior .* mejora_error_superior .*K_step+ 0.33* K_step .* .!mejora_error_superior) .* .!mejora_error_negativa + 0.33*K_step .*mejora_error_negativa
                #println(K_step)

                inc_valor_restricciones_no_fijas_total=inc_valor_restricciones_no_fijas_total-K_step.*gradiente_fuerzas_nodos

            end
            

            # if iter_gradiente==1
            # else
            #     for i in 1:length(inc_valor_restricciones_no_fijas_total_old)
            #         println(gdl_no_fijos_total[i]," ",i," ",gradiente_fuerzas_nodos[i]," ",gradiente_fuerzas_old[i]," ", mejora_error_inferior[i]," ", mejora_error_superior[i])
            #         println(gdl_no_fijos_total[i]," ",inc_valor_restricciones_no_fijas_total_old[i]," --> ",inc_valor_restricciones_no_fijas_total[i])
            #     end
            # end
            gradiente_fuerzas_old=gradiente_fuerzas_nodos
            println("error fuerzas acumulado entre todos los nodos= ", sum(abs.(gradiente_fuerzas_nodos)))
            push!(error_total, sum(abs.(gradiente_fuerzas_nodos)))
            push!(error_fuerzas,gradiente_fuerzas_nodos)
            display(plot(error_total, title="error acumulado entre todos los nodos"))
            display(plot(error_fuerzas, title="error en cada nodo"))
            # readline()
            

            ### RESUMEN DATOS ENERGIAS
            println()
            println("ENERGIA INCREMENTADA HEXAEDROS METAMATERIAL VS SOLIDO CONTINUO: ")
            for hexaedro in hexaedros_analizados_metam_nuevo 

                # si el elemento se estudia por primera vez en este incremento
                if (findall(in(hexaedro), hexaedros_analizados_metam_nuevo)!=[]) &&  (findall(in(hexaedro), hexaedros_analizados_metam)!=[])
                    println("ELEMENTO: ", hexaedro, " inc energia como continuo= ",inc_energias_solidos_gruesos[hexaedro], " inc energia como metam= ",inc_energia_elemento_solido_grueso_metam[hexaedro])
                else
                    println("ELEMENTO: ", hexaedro, " energia como continuo= ",energias_solidos_gruesos[hexaedro], " energia como metam= ",inc_energia_elemento_solido_grueso_metam[hexaedro])
                end
            end
            println()
            println("INC ENERGIA TOTAL SOLIDO CONTINUO ",sum(inc_energias_hexaedros))
            println("ENERGIA TOTAL SOLIDO CONTINUO ",sum(energias_hexaedros_acumulado_nuevo))

            # SI SE ANALIZAN TODOS LOS HEXAEDROS
            if length(hexaedros_analizados_metam_nuevo)==N_elementos_solidos_grueso
                println("INC ENERGIA TOTAL METAMATERIAL ",sum(values(inc_energia_elemento_solido_grueso_metam)))
                println("ENERGIA TOTAL METAMATERIAL ",sum(values(energia_elemento_solido_grueso_metam)))
            else
            end

            #readline() 
        end
    
       


        ### CRITERIO DE SALIDA
        SSR=0; N_hexaedros_analizados=length(hexaedros_analizados_metam_nuevo)
        for elemento_solido_grueso in hexaedros_analizados_metam_nuevo #collect(keys(inc_energia_elemento_solido_zona))
            # el criterio de salida será con sum of square residuals
            #SSR = SSR + 1/N_hexaedros_analizados*(inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]-inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso])^2
            SSR=SSR+abs(inc_energia_elemento_solido_grueso_metam[elemento_solido_grueso]-inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso])/inc_energias_hexaedros_grueso_aparente[elemento_solido_grueso] 
        end

        push!(SSR_vector, SSR)
        println();println("error SSR= ",SSR)
        if (SSR<tolerancia)#(abs(sqrt(SSR)-sqrt(SSR_antiguo))/sqrt(SSR) < tolerancia) || (abs(sqrt(SSR))/inc_energia_hexaedros_total<tolerancia_2)
            # Youngs_tangentes_hexaedros_finos=Youngs_tangentes_hexaedros_finos # este vector con los valores de los modulos de young lo usaremos para el siguiente valor de desplazamiento impuesto y asi facilitar la convergencia
            # println(Youngs_tangentes_hexaedros_finos)
            break
        else; 
            SSR_antiguo=SSR #;readline()
            if (iter==iter_max) && (abs(SSR) != minimum(abs.(SSR_vector)))
                println(SSR_vector)
                println("ATENCION: No ha habido buena convergencia: No se ha reducido el error entre el solido homogeneo y el metamaterial con esta ultima iteracion")
            else
            end
        end


        #readline()

        hexaedros_analizados_metam_nuevo==[]
    end


    return hexaedros_analizados_metam_nuevo, inc_u_solido, inc_f_solido, u_solido_acumulado_nuevo, f_solido_acumulado_nuevo, fuerzas_hexaedros_acumulado_nuevo, energias_hexaedros_acumulado_nuevo, inc_energia_interna_total_solido, inc_energias_solidos_gruesos, energias_solidos_gruesos, inc_deform_equiv_solido_grueso, inc_u_metam, inc_f_metam, f_acumulada_elementos_nuevo, Youngs_tangentes_hexaedros_finos, Youngs_tangentes_hexaedros_gruesos, Youngs_tangentes_memoria_metam_nuevo, epsilon_voigt_nuevo, energias_acumuladas_elementos_nuevo, eps_x_ptos_int_nuevo, eps_xy_ptos_int_nuevo, eps_xz_ptos_int_nuevo, gdl_restring_globales_metam, inc_valor_restricciones_fijos_global_metam, gdl_no_fijos_total, gdl_total_zona_hexaedros, gdl_fijos_zona_hexaedros, gdl_no_fijos_zona_hexaedros, gdl_total_hexaedro, gdl_fijos_hexaedro, gdl_restring_fijos_hexaedro_local, gdl_restring_no_fijos_hexaedro_local, gdl_contorno_solido_completo_no_fronteras, energia_elemento_solido_grueso_metam, inc_energia_elemento_solido_grueso_metam, invariante_def_memoria_metam_nuevo, inc_valor_restricciones_no_fijas_total, K_step, gradiente_fuerzas_nodos, datos_contribucion_barras_gdl
    
end #function resolucion_acoplamiento



end #module