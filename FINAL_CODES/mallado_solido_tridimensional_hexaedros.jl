#MODULO PARA CREAR LOS ELEMENTOS SOLIDOS Y GENERAR EL MALLADO (DISTINTAS POSIBILIDADES DE SOLIDOS: HEXAEDROS, TETRAEDROS, etc.)
module mallado_solido_tridimensional_hexaedros

using LinearAlgebra

export mallado_3D_hexaedros, volumenes_elementos_hexaedros, numeracion_caras_hexaedros ,hexaedros_zona_alrededor, determinar_caras_contorno_solido, is_in_hexaedro, is_in_alguna_cara_hexaedro, is_in_cara_hexaedro, tipo_caras_hexaedro, determinar_hexaedros_mallado_fino_en_grueso


function mallado_3D_hexaedros(X,Y,Z,Nx,Ny,Nz)

    size_x=X/Nx; Nodos_x=Nx+1
    size_y=Y/Ny; Nodos_y=Ny+1
    size_z=Z/Nz; Nodos_z=Nz+1
    N_nodos=Nodos_x*Nodos_y*Nodos_z
    N_elementos=Nx*Ny*Nz
    conectividad=zeros(Int64,N_elementos,9) # 4 nodos por elemento mas el numero de elemento
    
    # Definir el vector de coordenadas de los nodos: Se numeran los nodos recorriendo Y, luego X y luego Z
    coords=Array{Any}(zeros(N_nodos,4)) # 4 columnas: las 3 de coordenadas y el numero de nodo 
    n=1
    for k in 1:Nodos_z
        for i in 1:Nodos_x
            for j in 1:Nodos_y
                coords[n,1:3]=[(i-1)*size_x,(j-1)*size_y,(k-1)*size_z]
                coords[n,4]=n
                n=n+1
            end
        end
    end

    n=1
    # Determinar los nodos que corresponden a cada elemento
    for k in 1:Nodos_z-1
        for i in 1:Nodos_x-1
            for j in 1:Nodos_y-1
                # Primero obtenemos los 8 nodos que forman cada cubo en el que dividimos el solido (numeracion local, 3 indices)
                nodos_cubo_relativo=zeros(8,3); nodos_cubo_global=zeros(8)
                nodos_cubo_relativo[1,:]=[i,j,k]
                nodos_cubo_relativo[2,:]=[i,j+1,k]
                nodos_cubo_relativo[3,:]=[i+1,j,k]
                nodos_cubo_relativo[4,:]=[i+1,j+1,k]
                nodos_cubo_relativo[5,:]=[i,j,k+1]
                nodos_cubo_relativo[6,:]=[i,j+1,k+1]
                nodos_cubo_relativo[7,:]=[i+1,j,k+1]
                nodos_cubo_relativo[8,:]=[i+1,j+1,k+1]

                #Pasamos de numeracion local (3 indices) a la normal generica (solo un indice)
                for i in 1:8
                    nodos_cubo_global[i]=Nodos_x*Nodos_y*(nodos_cubo_relativo[i,3]-1)+Nodos_y*(nodos_cubo_relativo[i,1]-1)+nodos_cubo_relativo[i,2]
                end
                      
                conectividad[n,1:8]=nodos_cubo_global
                conectividad[n,9]=n
                n=n+1
                
            end
        end
    end
    conectividad=Array{Int64}(conectividad)

    return coords, conectividad

end #function mallado_3D_hexahedros


function hexaedros_zona_alrededor(conectividad_3D) # Determinar los hexaedros que rodean a otro hexaedro

    N_hexaedros=length(conectividad_3D[:,1])
    hexaedros_alrededor_1_capa=Dict()

    # 1 capa de hexaedros alrededor de otro
    for hexaedro in 1:N_hexaedros
        hexaedros_alrededor_1_capa[hexaedro]=[]

        for nodo_vertice in 1:8
            nodo=conectividad_3D[hexaedro, nodo_vertice]
            hexaedros_contacto=findall(in(nodo),conectividad_3D[:,1:8])

            for hexaedro_contacto in hexaedros_contacto
                push!(hexaedros_alrededor_1_capa[hexaedro], hexaedro_contacto[1])
            end
           
        end

        unique!(sort!(hexaedros_alrededor_1_capa[hexaedro]))
    end

    return hexaedros_alrededor_1_capa

end #function hexaedros_zona_alrededor



function volumenes_elementos_hexaedros(coords, conectividad) # obtener el volumenes de los hexaedros que forman la malla
    N_elementos=length(conectividad[:,1])
    volum_hexaedros=zeros(N_elementos)

    for n in 1:N_elementos
        # Coordenadas de los nodos que componenn un elemento
        coords_elemento=zeros(8,3) #4 nodos por elemento y 3 coordenadas
        for i in 1:8 
            nodo=conectividad[n,i]
            coords_elemento[i,:]=coords[nodo,1:3]
        end
        longitud_x=coords_elemento[3,1]-coords_elemento[1,1]
        longitud_y=coords_elemento[2,2]-coords_elemento[1,2]
        longitud_z=coords_elemento[5,3]-coords_elemento[1,3]

        volum_hexaedros[n]=longitud_x*longitud_y*longitud_z  # producto mixto para sacar el volumen de un hexaedro

    end

    return volum_hexaedros

end #volumenes_elementos_hexaedros


function is_in_hexaedro(coords_punto,coords_vertices_hexaedro,volumen_hexaedro) # para determinar si un punto esta dentro (incliyendo las caras) de un hexaedro o fuera
    # Para determinar si un nodo está fuera o dentro de un elemento vemos si el volumen de los hexaedros que forma ese elemento es mayo o igual al del elemento tetraedrico
    volum_descompuesto=0; caras_pertenece=[]            
    tolerancia_distancia=1e-12

    #obtener los volumenes que forma el nodo de metamaterial con los nodos del hexaedro (hay que sumar el volumen de los minihexaedros)
    auxiliar_nodos=[1 2 4 3; 1 5 6 2; 3 4 8 7; 1 3 7 5; 2 6 8 4; 5 7 8 6] # cada minivolumen lo calculamos con 3 puntos del hexaedro y el del nodo de metamaterial que analizamos
    for i in 1:6
        indice1=auxiliar_nodos[i,1]; indice2=auxiliar_nodos[i,2]; indice3=auxiliar_nodos[i,3]; indice4=auxiliar_nodos[i,4]
        vector_AB=coords_vertices_hexaedro[indice2,:]-coords_vertices_hexaedro[indice1,:]
        vector_AC=coords_vertices_hexaedro[indice4,:]-coords_vertices_hexaedro[indice1,:]
        vector_AP=coords_punto-coords_vertices_hexaedro[indice1,:]

        vector_normal_plano=cross(vector_AB,vector_AC)
        area_plano=norm(vector_normal_plano)
        vector_normal_plano=vector_normal_plano/area_plano

        distancia_punto_plano=abs(dot(vector_AP,vector_normal_plano))

        if distancia_punto_plano < tolerancia_distancia
            push!(caras_pertenece,i)    # para ver en que cara o caras(si esta en una arista) se situa el punto
        else
        end

        volum_descompuesto=volum_descompuesto+1/3*distancia_punto_plano*area_plano
    end

    tolerancia_volumen=1e-10
    if abs(volum_descompuesto-volumen_hexaedro)/volumen_hexaedro < tolerancia_volumen
        if caras_pertenece == []
            valor=1 #en este caso está dentro del hexaedro 
        else
            valor=2 #en este caso está en la cara del hexaedro
        end
    else
        valor=0 #en este caso esta fuera del volumen_hexaedro
        caras_pertenece =[]
    end
    #println(valor);println(caras_pertenece)

    return valor, caras_pertenece

end # function is_in_hexaedro


function numeracion_caras_hexaedros(conectividad_3D)
    auxiliar_caras=[1 2 4 3; 1 5 6 2; 3 4 8 7; 1 3 7 5; 2 6 8 4; 5 7 8 6]
    N_hexaedros=length(conectividad_3D[:,1])
    k=1

    caras_hexaedros=zeros(Int64,N_hexaedros,6); # aqui se almacenan las caras (con la numeracion global) que le corresponden a cada hexaedro: cada fila es un hexaedro
    vertices_caras=0

    for hexaedro in 1:N_hexaedros
        i=1   
        
        for cara in 1:6
            vertices_cara=conectividad_3D[hexaedro, auxiliar_caras[cara,:]]'
            #println(vertices_cara)

            if k>1
                # No hay que repetir las caras (con que 3 de los vertices esten en una cara, significa que los 4 estan y que por tanto, es repetida)
                indices_vertices_caras_1=findall(in(vertices_cara[1]),vertices_caras)
                
                # si encontramos al vertice 1 en alguna cara, seguimos
                if indices_vertices_caras_1!=[]
                    indices_vertices_caras_1=getindex.(indices_vertices_caras_1,1)
                    #println(indices_vertices_caras_1)
                    indices_vertices_caras_2=findall(in(vertices_cara[2]),vertices_caras[indices_vertices_caras_1,:])

                    # si encontramos al vertice 2 en alguna de las caras en la que esta el 1, seguimos
                    if indices_vertices_caras_2!=[]
                        indices_vertices_caras_2=getindex.(indices_vertices_caras_2,1)
                        #println(indices_vertices_caras_2)
                        indices_vertices_caras_3=findall(in(vertices_cara[3]),vertices_caras[indices_vertices_caras_1,:][indices_vertices_caras_2,:])

                        # si encontramos al vertice 3 en alguna de las caras en la que estan el 1 y el 2, con que 3 de los vertices esten en una cara, significa que los 4 estan y que por tanto, es repetida
                        if indices_vertices_caras_3!=[]
                            indices_vertices_caras_3=indices_vertices_caras_3[1][1]
                            cara_repetida=indices_vertices_caras_1[indices_vertices_caras_2[indices_vertices_caras_3]]
                            #println(cara_repetida)
                            caras_hexaedros[hexaedro,i]=cara_repetida
                            i=i+1
                        
                        # si el vertice 3 no esta en ninguna cara en la que tambien esten el 1 y el 2, significa que es nueva
                        else
                            caras_hexaedros[hexaedro,i]=k
                            vertices_caras=vcat(vertices_caras, vertices_cara)
                            k=k+1; i=i+1
                        end
                        
                    # si el vertice 2 no esta en ninguna cara en la que tambien esta el 1, significa que es nueva
                    else
                        caras_hexaedros[hexaedro,i]=k
                        vertices_caras=vcat(vertices_caras, vertices_cara)
                        k=k+1; i=i+1
                    end

                    # si el vertice 1 no esta en ninguna cara, significa que es nueva
                else
                    caras_hexaedros[hexaedro,i]=k
                    vertices_caras=vcat(vertices_caras, vertices_cara)
                    k=k+1; i=i+1
                end

            else
                caras_hexaedros[hexaedro,i]=k
                vertices_caras=vertices_cara
                k=k+1; i=i+1
            end
            
        end
    end

    return vertices_caras, caras_hexaedros

end # function


function determinar_caras_contorno_solido(X,Y,Z, vertices_caras, coords_3D) # Para ver si las caras del hexaedro estan en el exterior del solido o no
    
    N_caras=length(vertices_caras[:,1])
    epsilon=1e-12
    caras_contorno=[]

    for cara in 1:N_caras

        vertice_1=vertices_caras[cara,1]; vertice_2=vertices_caras[cara,2]; vertice_3=vertices_caras[cara,3];

        if abs(coords_3D[vertice_1,1]-0) < epsilon
            if abs(coords_3D[vertice_2,1]-0) < epsilon
                if abs(coords_3D[vertice_3,1]-0) < epsilon
                    push!(caras_contorno, cara); break
                else 
                end        
            else 
            end               
        else 
        end

        if abs(coords_3D[vertice_1,1]-X) < epsilon
            if abs(coords_3D[vertice_2,]-X) < epsilon
                if abs(coords_3D[vertice_3,1]-X) < epsilon
                    push!(caras_contorno, cara); break
                else 
                end        
            else 
            end               
        else 
        end

        if abs(coords_3D[vertice_1,2]-0) < epsilon
            if abs(coords_3D[vertice_2,2]-0) < epsilon
                if abs(coords_3D[vertice_3,2]-0) < epsilon
                    push!(caras_contorno, cara); break
                else 
                end        
            else 
            end               
        else 
        end

        if abs(coords_3D[vertice_1,2]-Y) < epsilon
            if abs(coords_3D[vertice_2,2]-Y) < epsilon
                if abs(coords_3D[vertice_3,2]-Y) < epsilon
                    push!(caras_contorno, cara); break
                else 
                end        
            else 
            end               
        else 
        end

        if abs(coords_3D[vertice_1,3]-0) < epsilon
            if abs(coords_3D[vertice_2,3]-0) < epsilon
                if abs(coords_3D[vertice_3,3]-0) < epsilon
                    push!(caras_contorno, cara); break
                else 
                end        
            else 
            end               
        else 
        end

        if abs(coords_3D[vertice_1,3]-Z) < epsilon
            if abs(coords_3D[vertice_2,3]-Z) < epsilon
                if abs(coords_3D[vertice_3,3]-Z) < epsilon
                    push!(caras_contorno, cara); break
                else 
                end        
            else 
            end               
        else 
        end

    end

    return caras_contorno
    
end


function tipo_caras_hexaedro(coords_vertices_hexaedro, X,Y,Z) # Para ver si las caras del hexaedro estan en el exterior del solido o no
    auxiliar_caras=[1 2 4 3; 1 5 6 2; 3 4 8 7; 1 3 7 5; 2 6 8 4; 5 7 8 6]
    epsilon=1e-12

    tipos_caras=zeros(Int64,6)

    # CARA 1 (BASE INFERIOR)
    if abs(coords_vertices_hexaedro[auxiliar_caras[1,1],3]-0) < epsilon
        tipos_caras[1]=1
    else 
    end

    # CARA 2 (CARA ANTERIOR)
    if abs(coords_vertices_hexaedro[auxiliar_caras[2,1],1]-0) < epsilon
        tipos_caras[2]=1
    else 
    end

    # CARA 3 (CARA POSTERIOR)
    if abs(coords_vertices_hexaedro[auxiliar_caras[3,1],1]-X) < epsilon
        tipos_caras[3]=1
    else 
    end

    # CARA 4 (LATERAL DERECHO visto desde atras)
    if abs(coords_vertices_hexaedro[auxiliar_caras[4,1],2]-0) < epsilon
        tipos_caras[4]=1
    else 
    end

    # CARA 5 (LATERAL IZQUIERDO visto desde atras)
    if abs(coords_vertices_hexaedro[auxiliar_caras[5,1],2]-Y) < epsilon
        tipos_caras[5]=1
    else 
    end
        
    # CARA 6 (BASE SUPERIOR)
    if abs(coords_vertices_hexaedro[auxiliar_caras[6,1],3]-Z) < epsilon
        tipos_caras[6]=1
    else 
    end

    return tipos_caras
    
end


function determinar_hexaedros_mallado_fino_en_grueso(N_3D_x_fino, N_3D_y_fino, N_3D_z_fino, N_3D_x_grueso, N_3D_y_grueso, N_3D_z_grueso)
    K_x=Int64(N_3D_x_fino/N_3D_x_grueso)
    K_y=Int64(N_3D_y_fino/N_3D_y_grueso)
    K_z=Int64(N_3D_z_fino/N_3D_z_grueso)

    N_hexaedros_gruesos=N_3D_x_grueso*N_3D_y_grueso*N_3D_z_grueso

    hexaedros_mallado_fino_en_grueso=zeros(Int64, N_hexaedros_gruesos, K_x*K_y*K_z)

    hexaedro_grueso=1
    for k in 1:N_3D_z_grueso
        for i in 1:N_3D_x_grueso
            for j in 1:N_3D_y_grueso              
                # obtenemos los hexaedros pequeños dentro de los grandes pero con una numeracion en coordenadas
                indices_x=collect((i-1)*K_x+1:i*K_x)
                indices_y=collect((j-1)*K_y+1:j*K_y)
                indices_z=collect((k-1)*K_z+1:k*K_z)

                m=1
                for indice_z in indices_z
                    for indice_x in indices_x
                        for indice_y in indices_y                        
                            hexaedros_mallado_fino_en_grueso[hexaedro_grueso, m]=indice_y + (indice_x-1)*N_3D_y_fino + (indice_z-1)*N_3D_y_fino*N_3D_x_fino
                            m=m+1
                        end
                    end
                end
                hexaedro_grueso=hexaedro_grueso+1
            end
        end
    end

    return hexaedros_mallado_fino_en_grueso

end


# function is_in_hexaedro(coords_punto,coords_vertices_hexaedro,volumen_hexaedro) # para determinar si un punto esta dentro (incliyendo las caras) de un hexaedro o fuera
#     # Para determinar si un nodo está fuera o dentro de un elemento vemos si el volumen de los hexaedros que forma ese elemento es mayo o igual al del elemento tetraedrico
#     volum_descompuesto=0            

#     #obtener los volumenes que forma el nodo de metamaterial con los nodos del hexaedro (hay que sumar el volumen de los minihexaedros)
#     auxiliar_nodos=[1 2 4 3; 1 5 6 2; 3 4 8 7; 1 3 7 5; 2 6 8 4; 5 7 8 6] # cada minivolumen lo calculamos con 3 puntos del hexaedro y el del nodo de metamaterial que analizamos
#     for i in 1:6
#         indice1=auxiliar_nodos[i,1]; indice2=auxiliar_nodos[i,2]; indice3=auxiliar_nodos[i,3]; indice4=auxiliar_nodos[i,4]
#         vector_AB=coords_vertices_hexaedro[indice2,:]-coords_vertices_hexaedro[indice1,:]
#         vector_AC=coords_vertices_hexaedro[indice4,:]-coords_vertices_hexaedro[indice1,:]
#         vector_AP=coords_punto-coords_vertices_hexaedro[indice1,:]

#         vector_normal_plano=cross(vector_AB,vector_AC)
#         area_plano=norm(vector_normal_plano)
#         vector_normal_plano=vector_normal_plano/area_plano

#         distancia_punto_plano=abs(dot(vector_AP,vector_normal_plano))

#         volum_descompuesto=volum_descompuesto+1/3*distancia_punto_plano*area_plano
#     end

#     tolerancia_volumen=1e-10
#     if abs(volum_descompuesto-volumen_hexaedro)/volumen_hexaedro < tolerancia_volumen
#         valor=true #en este caso está dentro del hexaedro (o en alguna cara)
#     else
#         valor=false #en este caso esta fuera del volumen_hexaedro
#     end

#     return valor

# end # function is_in_hexaedro


# function is_in_alguna_cara_hexaedro(coords_punto,coords_vertices_hexaedro) # funcion para ver si un punto esta en aguna de las caras del hexaedro
    
#     valor=false #esta variable determina si esta en alguna cara (true) o no (false)
#     tolerancia_area=1e-10         

#     #obtener los volumenes que forma el nodo de metamaterial con los nodos del hexaedro (hay que sumar el volumen de los minihexaedros)
#     auxiliar_nodos=[1 2 4 3; 1 5 6 2; 3 4 8 7; 1 3 7 5; 2 6 8 4; 5 7 8 6] # cada minivolumen lo calculamos con 3 puntos del hexaedro y el del nodo de metamaterial que analizamos
#     for i in 1:6
#         indice1=auxiliar_nodos[i,1]; indice2=auxiliar_nodos[i,2]; indice3=auxiliar_nodos[i,3]; indice4=auxiliar_nodos[i,4]
#         vector_AB=coords_vertices_hexaedro[indice2,:]-coords_vertices_hexaedro[indice1,:]
#         vector_AC=coords_vertices_hexaedro[indice4,:]-coords_vertices_hexaedro[indice1,:]
#         vector_AP=coords_punto-coords_vertices_hexaedro[indice1,:]

#         vector_normal_plano=cross(vector_AB,vector_AC)
#         area_plano=norm(vector_normal_plano)
#         vector_normal_plano=vector_normal_plano/area_plano

#         distancia_punto_plano=abs(dot(vector_AP,vector_normal_plano))

#         if distancia_punto_plano < tolerancia_area
#             valor=true #en este caso está dentro del hexaedro (o en alguna cara)
#             break
#         else
#         end
    
#     end

#     return valor

# end # function is_in_alguna_cara


# function is_in_cara_hexaedro(coords_punto,coords_vertices_cara_hexaedro) # funcion para ver si un punto esta en aguna de las caras del hexaedro
    
#     valor=false #esta variable determina si esta en alguna cara (true) o no (false)

#     vector_AB=coords_vertices_hexaedro[2,:]-coords_vertices_cara_hexaedro[1,:]
#     vector_AC=coords_vertices_hexaedro[4,:]-coords_vertices_cara_hexaedro[1,:]
#     vector_AP=coords_punto-coords_vertices_cara_hexaedro[1,:]

#     vector_normal_plano=cross(vector_AB,vector_AC)
#     area_plano=norm(vector_normal_plano)
#     vector_normal_plano=vector_normal_plano/area_plano

#     distancia_punto_plano=abs(dot(vector_AP,vector_normal_plano))

#     if distancia_punto_plano < tolerancia_area
#         valor=true #en este caso está dentro del hexaedro (o en alguna cara)
#     else
#     end
    
#     return valor

# end # function is_in_cara


# function tipo_caras_hexaedro(coords_vertices_hexaedro, X,Y,Z) # Para ver si las caras del hexaedro estan en el exterior del solido o no
#     auxiliar_caras=[1 2 4 3; 1 5 6 2; 3 4 8 7; 1 3 7 5; 2 6 8 4; 5 7 8 6]
#     epsilon=1e-12

#     tipos_caras=zeros(6)

#     # CARA 1 (BASE INFERIOR)
#     if abs(coords_vertices_hexaedro[auxiliar_caras[i,j],3]-0) < epsilon
#         tipos_caras[1]=1
#     else 
#     end

#     # CARA 2 (CARA ANTERIOR)
#     if abs(coords_vertices_hexaedro[auxiliar_caras[i,j],1]-0) < epsilon
#         tipos_caras[2]=1
#     else 
#     end

#     # CARA 3 (CARA POSTERIOR)
#     if abs(coords_vertices_hexaedro[auxiliar_caras[i,j],1]-X) < epsilon
#         tipos_caras[3]=1
#     else 
#     end

#     # CARA 4 (LATERAL DERECHO visto desde atras)
#     if abs(coords_vertices_hexaedro[auxiliar_caras[i,j],2]-0) < epsilon
#         tipos_caras[4]=1
#     else 
#     end

#     # CARA 5 (LATERAL IZQUIERDO visto desde atras)
#     if abs(coords_vertices_hexaedro[auxiliar_caras[i,j],2]-Y) < epsilon
#         tipos_caras[5]=1
#     else 
#     end
        
#     # CARA 6 (BASE SUPERIOR)
#     if abs(coords_vertices_hexaedro[auxiliar_caras[i,j],3]-Z) < epsilon
#         tipos_caras[6]=1
#     else 
#     end

#     return tipos_caras

# end



end # module