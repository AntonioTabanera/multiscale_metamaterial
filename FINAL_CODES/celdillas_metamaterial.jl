#MODULO PARA CREAR LAS CELDILLAS DE METAMATERIAL Y GENERAR EL MALLADO (DISTINTAS POSIBILIDADES DE CELDILLAS: FCC, BCC, etc.)
module celdillas_metamaterial    


export celdilla_FCC, celdilla_FCC_nodos_itermedios, gdl_metamaterial_contorno_solido_no_fronteras, celdilla_simple, factor_barras_caras_solido_gradiente

function celdilla_FCC(X,Y,Z,Nx,Ny,Nz)
    # X,Y,Z = Diensiones de la probeta
    # Nx,Ny,Nz = numero de celdillas en cada direccion

    x=X/Nx #dimension x celdilla
    y=Y/Ny #dimension x celdilla
    z=Z/Nz   #dimension x celdilla
    dim_celdilla=[x,y,z]

    #Capas_x=2*Nx+1; Capas_y=2*Ny+1; Capas_z=2*Nz+1;   # Capas en cada direccion 
    #Nodos_total=1/2*Capas_x*Capas_y*(Capas_z-1)+         #Nodos totales en la probeta
    
    ## COORDENADAS DE LOS NODOS:
    coords=Array{Any}(undef,0,3);                 #Matriz que contiene la coordenadas x,y,z de cada nodo y el numero del nodo
    Nodos_arista_XY=[]; Nodos_arista_XZ=[]; Nodos_arista_YZ=[]; Nodos_cara_X=[]; Nodos_cara_Y=[]; Nodos_cara_Z=[]; nodos_origen_celdilla=[]

    n=1; epsilon=1e-12
    # Se numeran los nodos recorriendo Y, luego X y luego Z
    for k in 1:Nz+1
        for i in 1:Nx+1
            for j in 1:Ny+1
                # global n  #(esto solo hace falta si en vez de en el módulo se pone directamente como parte del programa principal)
                coord_x=(i-1)*x; coord_y=(j-1)*y ; coord_z=(k-1)*z
                coords=vcat(coords,[coord_x  coord_y   coord_z])

                # clasificar al nodo en interior o perteneciente a alguna cara o arista
                if abs(coord_x-X) < epsilon 
                    if abs(coord_y-Y) < epsilon
                        if abs(coord_z-Z) < epsilon
                        else; push!(Nodos_arista_XY,n)
                        end
                    else 
                        if abs(coord_z-Z) < epsilon; push!(Nodos_arista_XZ,n)
                        else ; push!(Nodos_cara_X,n)
                        end   
                    end              
                else 
                    if abs(coord_y-Y) < epsilon
                        if abs(coord_z-Z) < epsilon; push!(Nodos_arista_YZ,n)
                        else; ; push!(Nodos_cara_Y,n)
                        end
                    else 
                        if abs(coord_z-Z) < epsilon; push!(Nodos_cara_Z,n)
                        else ; push!(nodos_origen_celdilla,n)
                        end                  
                    end
                end 
                n=n+1
            end
        end


    end
    N_nodos_total=n-1 #numero de nodos total de la estructura 
    coords=hcat(coords,1:N_nodos_total)

                    
    ## ESTABLECER LA conectividad
    conectividad=Array{Int64}(undef,0,3);
    
    n=1
    for i in nodos_origen_celdilla
        nodo_1=i;                nodo_2=i+1;            nodo_3=i+Ny+1;          nodo_4=nodo_3+1;     #nodos de la capa baja de la celdilla
        nodo_5=i+(Nx+1)*(Ny+1);  nodo_6=nodo_5+1;       nodo_7=nodo_5+Ny+1;     nodo_8=nodo_7+1;     #nodos de la capa superior de la celdilla

        
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_5 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_6 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_7 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_2 nodo_5 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_3 nodo_5 n]);    n=n+1
    
    end            
    for i in Nodos_cara_Y  #la cara que esta a una distancia Y del origen (vector normal el Y)
        nodo_1=i;               nodo_2=i+Ny+1;         nodo_3=i+(Nx+1)*(Ny+1);       nodo_4=nodo_3+Ny+1;        

        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1

    end
    for i in Nodos_cara_X  
        nodo_1=i;               nodo_2=i+1;         nodo_3=i+(Ny+1)*(Nx+1);       nodo_4=nodo_3+1 

        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1

    end
    for i in Nodos_cara_Z
        nodo_1=i;               nodo_2=i+1;         nodo_3=i+Ny+1;          nodo_4=nodo_3+1;      

        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1
        
    end
    for i in Nodos_arista_XY
        nodo_1=i;               nodo_2=i+(Nx+1)*(Ny+1)
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1

    end
    for i in Nodos_arista_XZ
        nodo_1=i;               nodo_2=i+1;
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1

    end
    for i in Nodos_arista_YZ
        nodo_1=i;               nodo_2=i+Ny+1;
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1

    end

    return coords, conectividad

end # celdilla_FCC


function celdilla_FCC_nodos_itermedios(X,Y,Z,Nx,Ny,Nz) #generar la geometria pero considerando los nodos en los que se cruzan las barras
    # X,Y,Z = Diensiones de la probeta
    # Nx,Ny,Nz = numero de celdillas en cada direccion

    x=X/Nx #dimension x celdilla
    y=Y/Ny #dimension x celdilla
    z=Z/Nz   #dimension x celdilla
    dim_celdilla=[x,y,z]

    #Capas_x=2*Nx+1; Capas_y=2*Ny+1; Capas_z=2*Nz+1;   # Capas en cada direccion 
    #Nodos_total=1/2*Capas_x*Capas_y*(Capas_z-1)+         #Nodos totales en la probeta
    
    ## COORDENADAS DE LOS NODOS:
    coords=Array{Any}(undef,0,3);                 #Matriz que contiene la coordenadas x,y,z de cada nodo y el numero del nodo
    Nodos_arista_XY=[]; Nodos_arista_XZ=[]; Nodos_arista_YZ=[]; Nodos_cara_X=[]; Nodos_cara_Y=[]; Nodos_cara_Z=[]; nodos_origen_celdilla=[]

    n=1; epsilon=1e-12
    # Se numeran los nodos recorriendo Y, luego X y luego Z
    for k in 1:2*Nz+1
        if k%2 != 0 #si es impar la capa en z será una de las bases de la celdilla
            for i in 1:2*Nx+1
                for j in 1:Ny+1/2*(1-(-1)^i)
                    # global n  #(esto solo hace falta si en vez de en el módulo se pone directamente como parte del programa principal)
                    coord_x=(i-1)*x/2; coord_y=(j-1)*y + 1/4*(1+ (-1)^i)*y; coord_z=(k-1)*z/2
                    coords=vcat(coords,[coord_x  coord_y   coord_z])

                    # clasificar al nodo en interior o perteneciente a alguna cara o arista
                    if i%2!=0 # esto se hace para que solo coincida con las caras de la celdilla
                        if abs(coord_x-X) < epsilon 
                            if abs(coord_y-Y) < epsilon
                                if abs(coord_z-Z) < epsilon
                                else; push!(Nodos_arista_XY,n)
                                end
                            else 
                                if abs(coord_z-Z) < epsilon; push!(Nodos_arista_XZ,n)
                                else ; push!(Nodos_cara_X,n)
                                end   
                            end              
                        else 
                            if abs(coord_y-Y) < epsilon
                                if abs(coord_z-Z) < epsilon; push!(Nodos_arista_YZ,n)
                                else; ; push!(Nodos_cara_Y,n)
                                end
                            else 
                                if abs(coord_z-Z) < epsilon; push!(Nodos_cara_Z,n)
                                else ; push!(nodos_origen_celdilla,n)
                                end                  
                            end
                        end 
                    else
                    end
                    n=n+1
                end
            end

        else #si es una capa par en Z, sus nodos no perteneceran a las caras superior o inferior de ninguna celdilla
            for i in 1:2*Nx+1
                for j in 1:Ny+1/2*(1+(-1)^i)
                    coord_x=(i-1)*x/2;  coord_y=(j-1)*y + 1/4*(1- (-1)^i)*y;  coord_z=(k-1)*z/2
                    coords=vcat(coords,[coord_x  coord_y   coord_z])
                    n=n+1
                end
            end
        end
    end
    N_nodos_total=n-1 #numero de nodos total de la estructura 
    coords=hcat(coords,1:N_nodos_total)

                    
    ## ESTABLECER LA conectividad
    conectividad=Array{Int64}(undef,0,3);
    
    n=1
    for i in nodos_origen_celdilla
        nodo_1=i;               nodo_2=i+1;         nodo_3=i+Ny+1;          nodo_4=nodo_3+Ny;        nodo_5=nodo_4+1 #nodos de la capa baja de la celdilla
        nodo_6=i+(Ny+1)*(Nx+1)+Ny*Nx;               nodo_7=nodo_6+Ny;       nodo_8=nodo_7+1;         nodo_9=nodo_8+Ny #nodos de la capa media de la celdilla
        nodo_10=nodo_6+(Ny)*(Nx+1)+(Ny+1)*Nx;       nodo_11=nodo_10+1;      nodo_12=nodo_10+Ny+1;    nodo_13=nodo_12+Ny;        nodo_14=nodo_13+1 #nodos de la capa alta de la celdilla        
        
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_6 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_7 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_10 n]);   n=n+1
        conectividad=vcat(conectividad,[nodo_3 nodo_5 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_6 nodo_11 n]);   n=n+1
        conectividad=vcat(conectividad,[nodo_7 nodo_13 n]);   n=n+1
        conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_2 nodo_6 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_3 nodo_4 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_6 nodo_10 n]);   n=n+1
        conectividad=vcat(conectividad,[nodo_4 nodo_7 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_7 nodo_10 n]);   n=n+1
    
    end            
    for i in Nodos_cara_Y  #la cara que esta a una distancia Y del origen (vector normal el Y)
        nodo_1=i;               nodo_2=i+2*Ny+1;         nodo_3=i+(Ny)*(Nx+1)+(Ny+1)*Nx+Ny+1;          
        nodo_4=nodo_3+Nx*Ny+Nx*(Ny+1);        nodo_5=nodo_4+2*Ny+1

        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_3 nodo_5 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_3 nodo_4 n]);    n=n+1

    end
    for i in Nodos_cara_X  
        nodo_1=i;               nodo_2=i+1;         nodo_3=i+(Ny+1)*(Nx+1)+Ny*Nx;          
        nodo_4=nodo_3+(Ny)*(Nx+1)+(Ny+1)*Nx;        nodo_5=nodo_4+1

        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_3 nodo_5 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_3 nodo_4 n]);    n=n+1

    end
    for i in Nodos_cara_Z
        nodo_1=i;               nodo_2=i+1;         nodo_3=i+Ny+1;          nodo_4=nodo_3+Ny;        nodo_5=nodo_4+1

        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_3 nodo_5 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_3 nodo_4 n]);    n=n+1
        
    end
    for i in Nodos_arista_XY
        nodo_1=i;               nodo_2=i+(Ny)*(Nx+1)+(Ny+1)*Nx+(Ny+1)*(Nx+1)+Ny*Nx;
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1

    end
    for i in Nodos_arista_XZ
        nodo_1=i;               nodo_2=i+1;
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1

    end
    for i in Nodos_arista_YZ
        nodo_1=i;               nodo_2=i+2*Ny+1;
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1

    end

    return coords, conectividad

end # celdilla_FCC_nodos_itermedios


function celdilla_simple(X,Y,Z,Nx,Ny,Nz)
    # X,Y,Z = Diensiones de la probeta
    # Nx,Ny,Nz = numero de celdillas en cada direccion

    x=X/Nx #dimension x celdilla
    y=Y/Ny #dimension x celdilla
    z=Z/Nz   #dimension x celdilla
    dim_celdilla=[x,y,z]

    #Capas_x=2*Nx+1; Capas_y=2*Ny+1; Capas_z=2*Nz+1;   # Capas en cada direccion 
    #Nodos_total=1/2*Capas_x*Capas_y*(Capas_z-1)+         #Nodos totales en la probeta
    
    ## COORDENADAS DE LOS NODOS:
    coords=Array{Any}(undef,0,3);                 #Matriz que contiene la coordenadas x,y,z de cada nodo y el numero del nodo
    Nodos_arista_XY=[]; Nodos_arista_XZ=[]; Nodos_arista_YZ=[]; Nodos_cara_X=[]; Nodos_cara_Y=[]; Nodos_cara_Z=[]; nodos_origen_celdilla=[]

    n=1; epsilon=1e-12
    # Se numeran los nodos recorriendo Y, luego X y luego Z
    for k in 1:Nz+1
        for i in 1:Nx+1
            for j in 1:Ny+1
                # global n  #(esto solo hace falta si en vez de en el módulo se pone directamente como parte del programa principal)
                coord_x=(i-1)*x; coord_y=(j-1)*y ; coord_z=(k-1)*z
                coords=vcat(coords,[coord_x  coord_y   coord_z])

                # clasificar al nodo en interior o perteneciente a alguna cara o arista
                if abs(coord_x-X) < epsilon 
                    if abs(coord_y-Y) < epsilon
                        if abs(coord_z-Z) < epsilon
                        else; push!(Nodos_arista_XY,n)
                        end
                    else 
                        if abs(coord_z-Z) < epsilon; push!(Nodos_arista_XZ,n)
                        else ; push!(Nodos_cara_X,n)
                        end   
                    end              
                else 
                    if abs(coord_y-Y) < epsilon
                        if abs(coord_z-Z) < epsilon; push!(Nodos_arista_YZ,n)
                        else; ; push!(Nodos_cara_Y,n)
                        end
                    else 
                        if abs(coord_z-Z) < epsilon; push!(Nodos_cara_Z,n)
                        else ; push!(nodos_origen_celdilla,n)
                        end                  
                    end
                end 
                n=n+1
            end
        end


    end
    N_nodos_total=n-1 #numero de nodos total de la estructura 
    coords=hcat(coords,1:N_nodos_total)

                    
    ## ESTABLECER LA conectividad
    conectividad=Array{Int64}(undef,0,3);
    
    n=1
    for i in nodos_origen_celdilla
        nodo_1=i;                nodo_2=i+1;            nodo_3=i+Ny+1;          nodo_4=nodo_3+1;     #nodos de la capa baja de la celdilla
        nodo_5=i+(Nx+1)*(Ny+1);  nodo_6=nodo_5+1;       nodo_7=nodo_5+Ny+1;     nodo_8=nodo_7+1;     #nodos de la capa superior de la celdilla

        
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_5 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_1 nodo_6 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_1 nodo_7 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_2 nodo_5 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_3 nodo_5 n]);    n=n+1
    
    end            
    for i in Nodos_cara_Y  #la cara que esta a una distancia Y del origen (vector normal el Y)
        nodo_1=i;               nodo_2=i+Ny+1;         nodo_3=i+(Nx+1)*(Ny+1);       nodo_4=nodo_3+Ny+1;        

        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1

    end
    for i in Nodos_cara_X  
        nodo_1=i;               nodo_2=i+1;         nodo_3=i+(Ny+1)*(Nx+1);       nodo_4=nodo_3+1 

        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1

    end
    for i in Nodos_cara_Z
        nodo_1=i;               nodo_2=i+1;         nodo_3=i+Ny+1;          nodo_4=nodo_3+1;      

        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1
        conectividad=vcat(conectividad,[nodo_1 nodo_3 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_1 nodo_4 n]);    n=n+1
        #conectividad=vcat(conectividad,[nodo_2 nodo_3 n]);    n=n+1
        
    end
    for i in Nodos_arista_XY
        nodo_1=i;               nodo_2=i+(Nx+1)*(Ny+1)
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1

    end
    for i in Nodos_arista_XZ
        nodo_1=i;               nodo_2=i+1;
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1

    end
    for i in Nodos_arista_YZ
        nodo_1=i;               nodo_2=i+Ny+1;
        conectividad=vcat(conectividad,[nodo_1 nodo_2 n]);    n=n+1

    end

    return coords, conectividad

end # celdilla_simple



# Determinar los nodos de metamaterial del contorno del solido 3D pero que no esten en aristas donde entren en contacto varios hexaedros
function gdl_metamaterial_contorno_solido_no_fronteras(X,Y,Z,N_3D_x,N_3D_y,N_3D_z,coords_metam)
    N_nodos=length(coords_metam[:,1])
    Nodos_contorno_solido_no_fronteras=[] #estan en el contorno del solido pero no en las fronteras de los hexaedros gruesos
    gdls_contorno_solido_no_fronteras=[]
    Nodos_contorno_solido_fronteras=[] #estan en el contorno del solido en las fronteras de los hexaedros gruesos
    gdls_contorno_solido_fronteras=[]

    # Determinar cual es el tamaño de los hexaedros
    inc_x=X/N_3D_x; inc_y=Y/N_3D_y; inc_z=Z/N_3D_z;

    epsilon=1e-13
    
    for nodo in 1:N_nodos
        coord_x=coords_metam[nodo,1]
        coord_y=coords_metam[nodo,2]
        coord_z=coords_metam[nodo,3]


        # Ver si el nodo esta en alguna cara y, de ser, en cual está
        esta_en_cara=[0,0,0]
        if abs(coord_x)<epsilon 
            esta_en_cara[1]=1
        else
        end

        if abs(coord_x-X)<epsilon
            esta_en_cara[1]=1
        else
        end

        if abs(coord_y)<epsilon 
            esta_en_cara[2]=1
        else
        end

        if abs(coord_y-Y)<epsilon
            esta_en_cara[2]=1
        else
        end

        if abs(coord_z)<epsilon 
            esta_en_cara[3]=1
        else
        end

        if abs(coord_z-Z)<epsilon
            esta_en_cara[3]=1
        else
        end
        

        # SI ESTA EN ALGUNA CARA
        # Ver si el nodo esta en alguna arista de los elementos en que se divide el solido y, de ser así, en cual está (siempre estará al menos en una arista por ser del contorno)
        if sum(esta_en_cara)!=0
            esta_en_arista=[0,0,0]
            if (coord_x % inc_x)<epsilon || inc_x-(coord_x % inc_x) <epsilon
                esta_en_arista[1]=1
            else 
            end

            if (coord_y % inc_y)<epsilon || inc_y-(coord_y % inc_y) <epsilon
                esta_en_arista[2]=1
            else 
            end

            if (coord_z % inc_z)<epsilon || inc_z-(coord_z % inc_z) <epsilon
                esta_en_arista[3]=1
            else 
            end
        else
        end
        #println(nodo);println(esta_en_cara);println(esta_en_arista);readline()

        
        # SI ESTA EN UNA SOLA CARA
        if sum(esta_en_cara)==1

            # SI ESTA SOLO EN UNA ARISTA (SIEMPRE ESTA AL MENOS EN UNA)
            if sum(esta_en_arista)==1
                append!(gdls_contorno_solido_no_fronteras,collect(6*nodo-5:6*nodo));push!(Nodos_contorno_solido_no_fronteras,nodo)

            # SI ESTA EN MAS DE UNA ARISTA
            else
                append!(gdls_contorno_solido_fronteras,collect(6*nodo-5:6*nodo));push!(Nodos_contorno_solido_fronteras,nodo)
            end

        # SI ESTA EN UNA ARISTA DEL SOLIDO COMPLETO (2 CARAS)
        elseif sum(esta_en_cara)==2

            # SI ESTA EN 3 ARISTAS A LA VEZ: SERÁ UNA ESQUINA DE ALGUN HEXAEDRO, NO SE CONSIDERA PARA EL CONTORNO
            if sum(esta_en_arista)==3          
                append!(gdls_contorno_solido_fronteras,collect(6*nodo-5:6*nodo));push!(Nodos_contorno_solido_fronteras,nodo)
            # SI NO ES UNA ESQUINA DE UN HEXAEDRO, SE CONSIDERA
            else
                append!(gdls_contorno_solido_no_fronteras,collect(6*nodo-5:6*nodo));push!(Nodos_contorno_solido_no_fronteras,nodo)
            end

        # SI ESTA EN UNA ESQUINA DEL SOLIDO COMPLETO (3 CARAS)
        elseif sum(esta_en_cara)==3
            append!(gdls_contorno_solido_no_fronteras,collect(6*nodo-5:6*nodo));push!(Nodos_contorno_solido_no_fronteras,nodo)

        # SI NO ESTA EN NINGUNA CARA, NO ES DEL CONTORNO
        else
        end



    end
    #println(Nodos_contorno_solido_no_fronteras)
    #println(length(Nodos_contorno_solido_no_fronteras))

    return Nodos_contorno_solido_no_fronteras, gdls_contorno_solido_no_fronteras, Nodos_contorno_solido_fronteras, gdls_contorno_solido_fronteras

end


function factor_barras_caras_solido_gradiente(X,Y,Z,conectividad_metam, coords_metam)

    factor_barras_caras_solido=ones(length(conectividad_metam[:,3]))

    epsilon=1e-13

    for i in conectividad_metam[:,3]
        nodo_cara=0; esta_en_cara=zeros(2,3)

        for k in 1:2

            nodo=conectividad_metam[i,k]

            coord_x=coords_metam[nodo,1]
            coord_y=coords_metam[nodo,2]
            coord_z=coords_metam[nodo,3]


            # Ver si el nodo esta en alguna cara y, de ser, en cual está
            
            if abs(coord_x)<epsilon 
                esta_en_cara[k,1]=1
            else
            end

            if abs(coord_x-X)<epsilon
                esta_en_cara[k,1]=1
            else
            end

            if abs(coord_y)<epsilon 
                esta_en_cara[k,2]=1
            else
            end

            if abs(coord_y-Y)<epsilon
                esta_en_cara[k,2]=1
            else
            end

            if abs(coord_z)<epsilon 
                esta_en_cara[k,3]=1
            else
            end

            if abs(coord_z-Z)<epsilon
                esta_en_cara[k,3]=1
            else
            end

        end

        if esta_en_cara!=zeros(2,3)
            indices_1=findall(in(1),esta_en_cara[1,:])
            indices_2=findall(in(1),esta_en_cara[2,:])

            caras_comunes=findall(in(indices_1), indices_2)

            if  length(caras_comunes)==1
                factor_barras_caras_solido[i]=2
            elseif length(caras_comunes)==2
                factor_barras_caras_solido[i]=4
            end
        else
        end
    end

    return factor_barras_caras_solido

end

end # module