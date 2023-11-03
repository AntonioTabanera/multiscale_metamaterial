#MODULO PARA REALIIZAR LOS CÁLCULOS MEF DE DEFORMACIONES, FUERZAS, DESPLAZAMIENTOS, ENERGÍA, etc. DEL SOLIDO (ELEMENTOS SOLIDOS TRIDIMENSIONALES)
module funciones_MEF_solido_hexaedros

include(raw"curva_sigma_epsilon_conjunto.jl"); using .curva_sigma_epsilon_conjunto
using LinearAlgebra
using StaticArrays
using SparseArrays
using DelimitedFiles

export  resolucion_solido_3D, condic_contorno_solido, desplazamientos_elemento, energias_y_fuerzas_acumuladas_elementos, resolucion_iterativa_solido


function matriz_K_elemento_hexaedro(coords_vertices_hexaedro, Young_elemento) # calcular la matriz de rigidez de un elemento
    poisson=Poisson_conjunto()
    K_elemento=zeros(24,24)
    chi_int=[-1/sqrt(3), 1/sqrt(3)]; eta_int=[-1/sqrt(3), 1/sqrt(3)]; tau_int=[-1/sqrt(3), 1/sqrt(3)];
    
    # Matriz constitutiva D: sigma=D*epsilon
    D=Young_elemento/((1+poisson)*(1-2*poisson))*[1-poisson poisson poisson 0 0 0;
                                                  poisson 1-poisson poisson 0 0 0;
                                                  poisson poisson 1-poisson 0 0 0;
                                                  0 0 0 (1-2*poisson)/2 0 0;
                                                  0 0 0 0 (1-2*poisson)/2 0;
                                                  0 0 0 0 0 (1-2*poisson)/2]

    for tau in tau_int
        for chi in chi_int
            for eta in eta_int
                # Vector de voigt de deformaciones con el siguente orden: [eps_x, eps_y, eps_z, gamma_xy, gamma_xz, gamma_yz]
                # DERIVADAS FUNCIONES DE FORMA RESPECTO A COORDENADAS ISOPARAMETRICAS: (coordenadas isoparametricas: chi, eta, tau)
                dN_dchi=zeros(8);                   dN_deta=zeros(8);                   dN_dtau=zeros(8)
                dN_dchi[1]=-1/8*(1-eta)*(1-tau);    dN_deta[1]=-1/8*(1-chi)*(1-tau);    dN_dtau[1]=-1/8*(1-chi)*(1-eta)
                dN_dchi[2]=-1/8*(1+eta)*(1-tau);    dN_deta[2]=1/8*(1-chi)*(1-tau);     dN_dtau[2]=-1/8*(1-chi)*(1+eta)
                dN_dchi[3]=1/8*(1-eta)*(1-tau);     dN_deta[3]=-1/8*(1+chi)*(1-tau);    dN_dtau[3]=-1/8*(1+chi)*(1-eta)
                dN_dchi[4]=1/8*(1+eta)*(1-tau);     dN_deta[4]=1/8*(1+chi)*(1-tau);     dN_dtau[4]=-1/8*(1+chi)*(1+eta)
                dN_dchi[5]=-1/8*(1-eta)*(1+tau);    dN_deta[5]=-1/8*(1-chi)*(1+tau);    dN_dtau[5]=1/8*(1-chi)*(1-eta)
                dN_dchi[6]=-1/8*(1+eta)*(1+tau);    dN_deta[6]=1/8*(1-chi)*(1+tau);     dN_dtau[6]=1/8*(1-chi)*(1+eta)
                dN_dchi[7]=1/8*(1-eta)*(1+tau);     dN_deta[7]=-1/8*(1+chi)*(1+tau);    dN_dtau[7]=1/8*(1+chi)*(1-eta)
                dN_dchi[8]=1/8*(1+eta)*(1+tau);     dN_deta[8]=1/8*(1+chi)*(1+tau);     dN_dtau[8]=1/8*(1+chi)*(1+eta)


                # MATRIZ JACOBIANA: (DERIVADAS DE COORDENADAS CARTESIANAS FRENTE A ISOPARAMETRICAS)
                dx_dchi=dot(dN_dchi,coords_vertices_hexaedro[:,1]);    dy_dchi=dot(dN_dchi,coords_vertices_hexaedro[:,2]);    dz_dchi=dot(dN_dchi,coords_vertices_hexaedro[:,3])
                dx_deta=dot(dN_deta,coords_vertices_hexaedro[:,1]);    dy_deta=dot(dN_deta,coords_vertices_hexaedro[:,2]);    dz_deta=dot(dN_deta,coords_vertices_hexaedro[:,3])
                dx_dtau=dot(dN_dtau,coords_vertices_hexaedro[:,1]);    dy_dtau=dot(dN_dtau,coords_vertices_hexaedro[:,2]);    dz_dtau=dot(dN_dtau,coords_vertices_hexaedro[:,3])
                
                Jacobiana=SMatrix{3, 3, Float64,9}([dx_dchi dx_deta dx_dtau;
                                                    dy_dchi dy_deta dy_dtau;
                                                    dz_dchi dz_deta dz_dtau])
                jacobiano=det(Jacobiana)


                # DERIVADAS FUNCIONES DE FORMA RESPECTO A COORDENADAS CARTESIANAS 
                inv_transp_Jac=inv(transpose(Jacobiana))
                dN_dx=zeros(8); dN_dy=zeros(8); dN_dz=zeros(8);
                for i in 1:8
                    dN_dx[i], dN_dy[i], dN_dz[i] = inv_transp_Jac * [dN_dchi[i], dN_deta[i], dN_dtau[i]]
                end


                # GENERAR MATRIZ B 
                B=zeros(6,24)
                for i in 1:8
                    B[:,3*i-2:3*i]=[dN_dx[i]    0        0;   
                                    0        dN_dy[i]    0;
                                    0           0     dN_dz[i]
                                    dN_dy[i] dN_dx[i]    0;   
                                    dN_dz[i]    0     dN_dx[i];
                                    0        dN_dz[i] dN_dy[i]]
                end

                # Matriz de rigidez del elemento K: F=K*u
                K_elemento=K_elemento+jacobiano*B'*D*B

            end
        end
    end

    return K_elemento

end # function matriz_K_elemento_tetraedro


function matriz_K_estructura_3D(coords_3D, conectividad_3D, Youngs_elementos)
    N_elementos=length(conectividad_3D[:,1]); N_nodos=length(coords_3D[:,1])
    gdl=3*N_nodos
    K_estructura=zeros(gdl,gdl)

    for elemento in 1:N_elementos
        coords_vertices_hexaedro=zeros(8,3);        gdl_elemento=[]

        for i in 1:8
            nodo=conectividad_3D[elemento,i]
            coords_vertices_hexaedro[i,:]=coords_3D[nodo,1:3]
            gdl_elemento=vcat(gdl_elemento,[3*nodo-2, 3*nodo-1, 3*nodo])
        end

        K_elemento=matriz_K_elemento_hexaedro(coords_vertices_hexaedro, Youngs_elementos[elemento])
        
        K_estructura[gdl_elemento,gdl_elemento]=K_estructura[gdl_elemento,gdl_elemento] + K_elemento
    end

    return K_estructura

end # matriz_K_estructura_3D


function condic_contorno_solido(punto_restriccion,gdl_restriccion,punto_desplazamiento,valor_desplazamiento,coords_3D) # condic de contorno en desplazamientos
    ## DETERMINA LOS GDL CON DESPLAZAMIENTO IMPUESTO Y CON RESTRICCIONES
    N_restricciones= length(punto_restriccion[:,1])
    N_desplazamientos= length(punto_desplazamiento[:,1])
    N_nodos= length(coords_3D[:,1])
    gdl_desplazamiento=valor_desplazamiento./valor_desplazamiento;  gdl_desplazamiento=replace(gdl_desplazamiento, NaN=>0)
    gdl_restring_globales=[];  gdl_desplaz_globales=[];  valor_desplaz_globales=[]#zeros(6*N_desplazamientos,2); 
    
    epsilon=1e-8
    coordenadas_x=sort(unique(coords_3D[:,1])); coordenadas_y=sort(unique(coords_3D[:,2])); coordenadas_z=sort(unique(coords_3D[:,3]))
    dimension_x_minima=coordenadas_x[2]-coordenadas_x[1];
    dimension_y_minima=coordenadas_y[2]-coordenadas_y[1];
    dimension_z_minima=coordenadas_z[2]-coordenadas_z[1];

    # RESTRICCIONES
    # se recorren todos los nodos y se comprueba se se situan en una coordenada donde haya impuesta una condicion de contorno
    for n in 1:N_restricciones

        a=length(gdl_restriccion) # esto se hace para ver si la condicion de contorno no coincide exactamente con ningun nodo
        for m in 1:N_nodos 

            gdlx=Int64(gdl_restriccion[n,1]*(3*m-2))
            gdly=Int64(gdl_restriccion[n,2]*(3*m-1))
            gdlz=Int64(gdl_restriccion[n,3]*(3*m))
            gdl_restring_globales_nodo=[gdlx,gdly,gdlz]

            if punto_restriccion[n,1]=="all" || abs(coords_3D[m,1]-punto_restriccion[n,1]) < epsilon
                if punto_restriccion[n,2]=="all"  || abs(coords_3D[m,2]-punto_restriccion[n,2]) < epsilon
                    if punto_restriccion[n,3]=="all" || abs(coords_3D[m,3]-punto_restriccion[n,3]) < epsilon
                        append!(gdl_restring_globales,gdl_restring_globales_nodo)
                        #println(coords_3D[m,:])
                    end
                else
                end
            else
            end

        end

        b=length(gdl_restring_globales) # si el vector de gdl de desplazamiento no ha crecido, significa que no se ha encontrado un nodo exacto
        
        ### SI NO SE HA ENCONTRADO NODO, SE LE IMPONE LA CONDICION A LOS MAS CERCANOS
        if b==a # si no se ha encontrado desplazamiento
            println("no se ha encontrado nodo, se impondrá la restriccion al mas cercano")
            # Esta aproximacion al nodo mas cercano será mas cierta cuanto mas celdillas haya
            for m in 1:N_nodos 

                gdlx=Int64(gdl_restriccion[n,1]*(3*m-2))
                gdly=Int64(gdl_restriccion[n,2]*(3*m-1))
                gdlz=Int64(gdl_restriccion[n,3]*(3*m))
                gdl_restring_globales_nodo=[gdlx,gdly,gdlz]
    
                if punto_restriccion[n,1]=="all" || abs(coords_3D[m,1]-punto_restriccion[n,1]) < dimension_x_minima
                    if punto_restriccion[n,2]=="all"  || abs(coords_3D[m,2]-punto_restriccion[n,2]) < dimension_y_minima
                        if punto_restriccion[n,3]=="all" || abs(coords_3D[m,3]-punto_restriccion[n,3]) < dimension_z_minima
                            append!(gdl_restring_globales,gdl_restring_globales_nodo)
                            #println(coords_3D[m,:])
                        end
                    else
                    end
                else
                end
    
            end
        else
        end
    end

    # DESPLAZAMIENTOS
    for n in 1:N_desplazamientos
        a=length(gdl_desplaz_globales) # esto se hace para ver si la condicion de contorno no coincide exactamente con ningun nodo

        for m in 1:N_nodos 

            gdlx=Int64(gdl_desplazamiento[n,1]*(3*m-2));           valor_desplaz_x=valor_desplazamiento[n,1]
            gdly=Int64(gdl_desplazamiento[n,2]*(3*m-1));           valor_desplaz_y=valor_desplazamiento[n,2]
            gdlz=Int64(gdl_desplazamiento[n,3]*(3*m));           valor_desplaz_z=valor_desplazamiento[n,3]
            gdl_desplaz_nodo=[gdlx,gdly,gdlz]
            valor_desplaz_nodo=[valor_desplaz_x,valor_desplaz_y,valor_desplaz_z]

            if punto_desplazamiento[n,1]=="all" || abs(coords_3D[m,1]-punto_desplazamiento[n,1]) < epsilon
                if punto_desplazamiento[n,2]=="all"  || abs(coords_3D[m,2]-punto_desplazamiento[n,2]) < epsilon
                    if punto_desplazamiento[n,3]=="all" || abs(coords_3D[m,3]-punto_desplazamiento[n,3]) < epsilon
                        append!(gdl_desplaz_globales, gdl_desplaz_nodo)
                        append!(valor_desplaz_globales, valor_desplaz_nodo)
                        #i=i+1
                        #println(coords_3D[m,:])
                    else
                    end
                else
                end
            else
            end

        end

        b=length(gdl_desplaz_globales) # si el vector de gdl de desplazamiento no ha crecido, significa que no se ha encontrado un nodo exacto
        
        ### SI NO SE HA ENCONTRADO NODO, SE LE IMPONE LA CONDICION A LOS MAS CERCANOS
        if b==a # si no se ha encontrado desplazamiento
            println("no se ha encontrado nodo, se impondrá el desplazamiento al mas cercano")
            # Esta aproximacion al nodo mas cercano será mas cierta cuanto mas divisiones haya
            for m in 1:N_nodos 

                gdlx=Int64(gdl_desplazamiento[n,1]*(3*m-2));           valor_desplaz_x=valor_desplazamiento[n,1]
                gdly=Int64(gdl_desplazamiento[n,2]*(3*m-1));           valor_desplaz_y=valor_desplazamiento[n,2]
                gdlz=Int64(gdl_desplazamiento[n,3]*(3*m));           valor_desplaz_z=valor_desplazamiento[n,3]
                gdl_desplaz_nodo=[gdlx,gdly,gdlz]
                valor_desplaz_nodo=[valor_desplaz_x,valor_desplaz_y,valor_desplaz_z]
    
                if punto_desplazamiento[n,1]=="all" || abs(coords_3D[m,1]-punto_desplazamiento[n,1]) < dimension_x_minima
                    if punto_desplazamiento[n,2]=="all"  || abs(coords_3D[m,2]-punto_desplazamiento[n,2]) < dimension_y_minima
                        if punto_desplazamiento[n,3]=="all" || abs(coords_3D[m,3]-punto_desplazamiento[n,3]) < dimension_z_minima
                            append!(gdl_desplaz_globales, gdl_desplaz_nodo)
                            append!(valor_desplaz_globales, valor_desplaz_nodo)
                            #i=i+1
                            #println(coords_3D[m,:])
                        else
                        end
                    else
                    end
                else
                end
    
            end
        else
        end

    end
    

    #Hay que eliminar los gdl de los nodos a los que se les ha impuesto una condicion porque no en todos los gdl tienen una
    gdl_no_nulo=findall(!in(0),gdl_restring_globales); gdl_restring_globales=gdl_restring_globales[gdl_no_nulo]
    gdl_no_nulo=findall(!in(0),valor_desplaz_globales); 
    gdl_desplaz_globales=gdl_desplaz_globales[gdl_no_nulo]
    valor_desplaz_globales=valor_desplaz_globales[gdl_no_nulo]
    

    if length(gdl_restring_globales)==0 && N_restricciones!=0
        println("no se ha podido encontrar ni aproximar ningun nodo para la restriccion")
    else
    end

    if length(gdl_desplaz_globales)==0 && N_desplazamientos!=0
        println("no se ha podido encontrar ni aproximar ningun nodo para el desplazamiento")
    else
    end


    # OBTENER TODOS LOS GDL CON ALGUNA RESTRICCION
    N_restricciones_movimiento_impedido=length(gdl_restring_globales)
    append!(gdl_restring_globales, gdl_desplaz_globales) #todos los nodos con ALGUNA condicion de contorno (y ordenados)
    valor_restricciones_globales=append!(zeros(N_restricciones_movimiento_impedido), valor_desplaz_globales)


    # Ver si a algun grado de libertad se repite porque se le han aplicado 2 condiciones de contorno diferentes
    if length(gdl_restring_globales)>length(unique(gdl_restring_globales))
        println("2 condiciones distintas para un mismo gdl")
    else
    end

    return gdl_restring_globales, valor_restricciones_globales

end # function condic_contorno_solido


function solver_desplazamientos_solido(gdl_restring_globales, valor_restricciones_globales, K_global) #obtener el vector u de desplazamientos globales
    gdl=length(K_global[:,1]); N_nodos=gdl/3
    gdl_totales=collect(1:gdl)
    gdl_libres=findall(!in(gdl_restring_globales),gdl_totales) #gdl sin restricciones
    u=zeros(gdl); f=zeros(gdl)
    u[gdl_restring_globales]=valor_restricciones_globales
    #u[gdl_libres]=K_global[gdl_libres,gdl_libres]\(-K_global[gdl_libres,gdl_restring_globales]*u[gdl_restring_globales])

    dimension_problema=length(gdl_libres);println("la dimension del problema del solido es: ",dimension_problema);
    dimension_limite=1000
    iter_max=10000

    if dimension_problema <= dimension_limite # con dimensiones pequeñas se resuelve la inversa numericamente

        u[gdl_libres]=K_global[gdl_libres,gdl_libres]\(-K_global[gdl_libres,gdl_restring_globales]*u[gdl_restring_globales]) #desplazamiento de los puntos sin ninguna restriccion
    
    else # con dimensiones grandes se resuelve la inversa iterativamente

        A=sparse(K_global[gdl_libres,gdl_libres]); b=-K_global[gdl_libres,gdl_restring_globales]*u[gdl_restring_globales]

        x=zeros(dimension_problema)
        residuo_antiguo=b; direccion=residuo_antiguo; k=0; norma_antigua=0; tolerancia=1e-10
        
        for i in 1:iter_max    
            alpha=(residuo_antiguo'*residuo_antiguo)/(direccion'*A*direccion)
            x=x+alpha*direccion
            norma_nueva=norm(x)

            # criterio de convergencia
            if abs(norma_antigua-norma_nueva)/norma_nueva < tolerancia;  break
            else
                residuo=residuo_antiguo-alpha*A*direccion
                beta=(residuo'*residuo)/(residuo_antiguo'*residuo_antiguo)
                direccion=residuo+beta*direccion
                residuo_antiguo=residuo
                norma_antigua=norma_nueva
            end
        end

        u[gdl_libres]=x

    end
    f[gdl_restring_globales]=K_global[gdl_restring_globales,gdl_libres]*u[gdl_libres]+K_global[gdl_restring_globales,gdl_restring_globales]*u[gdl_restring_globales]
    
    
    return u,f

end # function solver_desplazamientos_solido


function desplazamientos_elemento(coords_punto, elemento_solido, coords_3D, conectividad_3D, u_solido) #calcular interpolando los desplazamientos en un punto de un elemento
  
    x=coords_punto[1]; y=coords_punto[2]; z=coords_punto[3]
    u_punto=zeros(3) # los desplazamientos del punto

    u_elemento=zeros(8); v_elemento=zeros(8); w_elemento=zeros(8); # Vectores desplazamientos en los nodos del elemento
    for i in 1:8 
        nodo=conectividad_3D[elemento_solido,i]
        u_elemento[i]=u_solido[3*nodo-2]
        v_elemento[i]=u_solido[3*nodo-1]
        w_elemento[i]=u_solido[3*nodo]
    end

    #println("desplazamientos tetraedro =",u_elemento)
    vertice_1=conectividad_3D[elemento_solido,1]
    vertice_2=conectividad_3D[elemento_solido,2]
    vertice_3=conectividad_3D[elemento_solido,3]
    vertice_5=conectividad_3D[elemento_solido,5]

    longitud_x=coords_3D[vertice_3,1]-coords_3D[vertice_1,1]
    longitud_y=coords_3D[vertice_2,2]-coords_3D[vertice_1,2]
    longitud_z=coords_3D[vertice_5,3]-coords_3D[vertice_1,3]
    
    chi=-1+2*(x-coords_3D[vertice_1,1])/longitud_x
    eta=-1+2*(y-coords_3D[vertice_1,2])/longitud_y
    tau=-1+2*(z-coords_3D[vertice_1,3])/longitud_z

    # VALOR DE LAS FUNCIONES DE FORMA EN ESE PUNTO
    N_punto=zeros(8)
    N_punto[1]=1/8*(1-chi)*(1-eta)*(1-tau)    
    N_punto[2]=1/8*(1-chi)*(1+eta)*(1-tau)    
    N_punto[3]=1/8*(1+chi)*(1-eta)*(1-tau)    
    N_punto[4]=1/8*(1+chi)*(1+eta)*(1-tau)    
    N_punto[5]=1/8*(1-chi)*(1-eta)*(1+tau)  
    N_punto[6]=1/8*(1-chi)*(1+eta)*(1+tau) 
    N_punto[7]=1/8*(1+chi)*(1-eta)*(1+tau)     
    N_punto[8]=1/8*(1+chi)*(1+eta)*(1+tau)


    # Desplazamientos interpolados en un elemento
    u_punto[1]=N_punto'*u_elemento
    u_punto[2]=N_punto'*v_elemento 
    u_punto[3]=N_punto'*w_elemento

    return u_punto

end # function desplazamientos_elemento


function invariantes_deformaciones_elementos_hexaedricos(tipo_invariante, coords_3D, conectividad_3D, inc_u, epsilon_voigt_acumulado, Youngs_elementos) # Calculo en todos los elementos; al ser un elemento tetraedrico, solo hay un punto de integracion y solo un valor de las tensiones en todo el elemento
    poisson=Poisson_conjunto(); N_elementos=length(conectividad_3D[:,1])
    epsilon_voigt_nuevo=zeros(6,8,N_elementos)
    invariante_ptos_int=zeros(8) ; inc_def_equivalente_ptos_int=zeros(8) #; invariante_2=zeros(N_elementos); invariante_3=zeros(N_elementos); deformacion_equivalente=zeros(N_elementos)
    invariante=zeros(N_elementos) ; inc_def_equivalente=zeros(N_elementos)

    chi_int=[-1/sqrt(3), 1/sqrt(3)]; eta_int=[-1/sqrt(3), 1/sqrt(3)]; tau_int=[-1/sqrt(3), 1/sqrt(3)];
    

    for elemento in 1:N_elementos # recorremos todos los puntos de integracion de cada hexahedro

        coords_vertices_hexaedro=zeros(8,3); inc_u_elemento=zeros(24) #Vector desplazamientos en los nodos del elemento
        for i in 1:8 
            nodo=conectividad_3D[elemento,i]
            inc_u_elemento[3*i-2:3*i]=inc_u[3*nodo-2:3*nodo]
            coords_vertices_hexaedro[i,:]=coords_3D[nodo,1:3]
        end

        # Matriz constitutiva D: sigma=D*epsilon
        D=Youngs_elementos[elemento]/((1+poisson)*(1-2*poisson))*[1-poisson poisson poisson 0 0 0;
                                                        poisson 1-poisson poisson 0 0 0;
                                                        poisson poisson 1-poisson 0 0 0;
                                                        0 0 0 (1-2*poisson)/2 0 0;
                                                        0 0 0 0 (1-2*poisson)/2 0;
                                                        0 0 0 0 0 (1-2*poisson)/2]

        pto_int=0
        for tau in tau_int
            for chi in chi_int
                for eta in eta_int

                    pto_int=pto_int+1

                    # Vector de voigt de deformaciones con el siguente orden: [eps_x, eps_y, eps_z, gamma_xy, gamma_xz, gamma_yz]
                    # DERIVADAS FUNCIONES DE FORMA RESPECTO A COORDENADAS ISOPARAMETRICAS: (coordenadas isoparametricas: chi, eta, tau)
                    dN_dchi=zeros(8);                   dN_deta=zeros(8);                   dN_dtau=zeros(8)
                    dN_dchi[1]=-1/8*(1-eta)*(1-tau);    dN_deta[1]=-1/8*(1-chi)*(1-tau);    dN_dtau[1]=-1/8*(1-chi)*(1-eta)
                    dN_dchi[2]=-1/8*(1+eta)*(1-tau);    dN_deta[2]=1/8*(1-chi)*(1-tau);     dN_dtau[2]=-1/8*(1-chi)*(1+eta)
                    dN_dchi[3]=1/8*(1-eta)*(1-tau);     dN_deta[3]=-1/8*(1+chi)*(1-tau);    dN_dtau[3]=-1/8*(1+chi)*(1-eta)
                    dN_dchi[4]=1/8*(1+eta)*(1-tau);     dN_deta[4]=1/8*(1+chi)*(1-tau);     dN_dtau[4]=-1/8*(1+chi)*(1+eta)
                    dN_dchi[5]=-1/8*(1-eta)*(1+tau);    dN_deta[5]=-1/8*(1-chi)*(1+tau);    dN_dtau[5]=1/8*(1-chi)*(1-eta)
                    dN_dchi[6]=-1/8*(1+eta)*(1+tau);    dN_deta[6]=1/8*(1-chi)*(1+tau);     dN_dtau[6]=1/8*(1-chi)*(1+eta)
                    dN_dchi[7]=1/8*(1-eta)*(1+tau);     dN_deta[7]=-1/8*(1+chi)*(1+tau);    dN_dtau[7]=1/8*(1+chi)*(1-eta)
                    dN_dchi[8]=1/8*(1+eta)*(1+tau);     dN_deta[8]=1/8*(1+chi)*(1+tau);     dN_dtau[8]=1/8*(1+chi)*(1+eta)


                    # MATRIZ JACOBIANA: (DERIVADAS DE COORDENADAS CARTESIANAS FRENTE A ISOPARAMETRICAS)
                    dx_dchi=dot(dN_dchi,coords_vertices_hexaedro[:,1]);    dy_dchi=dot(dN_dchi,coords_vertices_hexaedro[:,2]);    dz_dchi=dot(dN_dchi,coords_vertices_hexaedro[:,3])
                    dx_deta=dot(dN_deta,coords_vertices_hexaedro[:,1]);    dy_deta=dot(dN_deta,coords_vertices_hexaedro[:,2]);    dz_deta=dot(dN_deta,coords_vertices_hexaedro[:,3])
                    dx_dtau=dot(dN_dtau,coords_vertices_hexaedro[:,1]);    dy_dtau=dot(dN_dtau,coords_vertices_hexaedro[:,2]);    dz_dtau=dot(dN_dtau,coords_vertices_hexaedro[:,3])
                         
                    Jacobiana=SMatrix{3, 3, Float64,9}([dx_dchi dx_deta dx_dtau;
                                                        dy_dchi dy_deta dy_dtau;
                                                        dz_dchi dz_deta dz_dtau])
                    jacobiano=det(Jacobiana)


                    # DERIVADAS FUNCIONES DE FORMA RESPECTO A COORDENADAS CARTESIANAS 
                    inv_transp_Jac=inv(transpose(Jacobiana))
                    dN_dx=zeros(8); dN_dy=zeros(8); dN_dz=zeros(8);
                    for i in 1:8
                        dN_dx[i], dN_dy[i], dN_dz[i] = inv_transp_Jac * [dN_dchi[i], dN_deta[i], dN_dtau[i]]
                    end


                    # GENERAR MATRIZ B 
                    B=zeros(6,24)
                    for i in 1:8
                        B[:,3*i-2:3*i]=[dN_dx[i]    0        0;   
                                        0        dN_dy[i]    0;
                                        0           0     dN_dz[i]
                                        dN_dy[i] dN_dx[i]    0;   
                                        dN_dz[i]    0     dN_dx[i];
                                        0        dN_dz[i] dN_dy[i]]
                    end


                    #Pasar de incrementos de deformaciones a deformaciones completas en cada punto de integracion
                    inc_epsilon_voigt=B*inc_u_elemento
                    epsilon_voigt = epsilon_voigt_acumulado[:,pto_int,elemento] + inc_epsilon_voigt
                    epsilon_voigt_nuevo[:,pto_int,elemento]=epsilon_voigt

                    # Tensor de deformaciones
                    epsilon_tensor=SMatrix{3, 3, Float64,9}([epsilon_voigt[1] epsilon_voigt[4]/2 epsilon_voigt[6]/2;   
                                                            epsilon_voigt[4]/2 epsilon_voigt[2] epsilon_voigt[5]/2;
                                                            epsilon_voigt[6]/2 epsilon_voigt[5]/2 epsilon_voigt[3]])

                    # Invariantes del tensor
                    inv_1=epsilon_voigt[1]+epsilon_voigt[2]+epsilon_voigt[3]       
                    signo=sign(inv_1)
                    inv_2=epsilon_voigt[1]*epsilon_voigt[2]+epsilon_voigt[1]*epsilon_voigt[3]+epsilon_voigt[2]*epsilon_voigt[3]-1/4*(epsilon_voigt[4]^2+epsilon_voigt[5]^2+epsilon_voigt[6]^2)       
                    inv_3=det(epsilon_tensor)
                    def_equiv=sqrt(1/((1+poisson)*(1-2*poisson))*((1-poisson)*(epsilon_voigt[1]^2+epsilon_voigt[2]^2+epsilon_voigt[3]^2)+2*poisson*(epsilon_voigt[1]*epsilon_voigt[2]+epsilon_voigt[1]*epsilon_voigt[3]+epsilon_voigt[2]*epsilon_voigt[3]))+1/(2*(1+poisson))*(epsilon_voigt[4]^2+epsilon_voigt[5]^2+epsilon_voigt[6]^2))
                    def_equiv_2=1/(sqrt(2)*(1+poisson))*sqrt((epsilon_voigt[1]-epsilon_voigt[2])^2+(epsilon_voigt[2]-epsilon_voigt[3])^2+(epsilon_voigt[3]-epsilon_voigt[1])^2+3/2*(epsilon_voigt[4]^2+epsilon_voigt[5]^2+epsilon_voigt[6]^2))


                    if tipo_invariante == "inv_1"
                    invariante_ptos_int[pto_int]=inv_1
                    elseif tipo_invariante == "inv_2"
                    invariante_ptos_int[pto_int]=signo*sqrt(abs(inv_2))
                    elseif tipo_invariante == "inv_3"
                    invariante_ptos_int[pto_int]=signo*cbrt(abs(inv_3))
                    elseif tipo_invariante == "deform_equiv"
                    invariante_ptos_int[pto_int]=signo*def_equiv
                    elseif tipo_invariante == "deform_equiv_2"
                    invariante_ptos_int[pto_int]=signo*def_equiv_2
                    else
                    println("No se ha introducido bien el tipo de invariante")
                    end

                    # Incrementos de invariantes en esta iteracion
                    inc_def_equivalente_ptos_int[pto_int]=sqrt(1/((1+poisson)*(1-2*poisson))*((1-poisson)*(inc_epsilon_voigt[1]^2+inc_epsilon_voigt[2]^2+inc_epsilon_voigt[3]^2)+2*poisson*(inc_epsilon_voigt[1]*inc_epsilon_voigt[2]+inc_epsilon_voigt[1]*inc_epsilon_voigt[3]+inc_epsilon_voigt[2]*inc_epsilon_voigt[3]))+1/(2*(1+poisson))*(inc_epsilon_voigt[4]^2+inc_epsilon_voigt[5]^2+inc_epsilon_voigt[6]^2))
                    #inc_def_equivalente[elemento]=1/(sqrt(2)*(1+Poisson_conjunto()))*sqrt((inc_epsilon_voigt[1]-inc_epsilon_voigt[2])^2+(inc_epsilon_voigt[2]-inc_epsilon_voigt[3])^2+(inc_epsilon_voigt[3]-inc_epsilon_voigt[1])^2+3/2*(inc_epsilon_voigt[4]^2+inc_epsilon_voigt[5]^2+inc_epsilon_voigt[6]^2))

                end
            end
        end
        invariante[elemento]=sqrt(sum(invariante_ptos_int.^2)/8);  # cogemos como valor de deformacion equivalente a la medis de los ptos de integracion
        inc_def_equivalente[elemento]=sqrt(sum(inc_def_equivalente_ptos_int.^2)/8)#sum(inc_def_equivalente_ptos_int)/8;#
        #println(invariante_ptos_int)
    end
    #println(invariante_ptos_int)
    return epsilon_voigt_nuevo, invariante, inc_def_equivalente

end # function invariantes_deformaciones_elemento


function inc_energias_elementos(inc_u, coords_3D, conectividad_3D, Youngs_elementos)
    N_elementos=length(conectividad_3D[:,1])
    inc_energ_elementos=zeros(N_elementos)
    

    for elemento in 1:N_elementos
        inc_u_elemento=zeros(24); coords_vertices_hexaedro=zeros(8,3) # Vector desplazamientos en los nodos del elemento
        for i in 1:8
            nodo=conectividad_3D[elemento,i]
            coords_vertices_hexaedro[i,:]=coords_3D[nodo,1:3]
            inc_u_elemento[3*i-2:3*i]=inc_u[3*nodo-2:3*nodo]
        end

        K_elemento=matriz_K_elemento_hexaedro(coords_vertices_hexaedro, Youngs_elementos[elemento])

        inc_energ_elementos[elemento]=1/2*inc_u_elemento'*K_elemento*inc_u_elemento

    end

    return inc_energ_elementos

end # function fuerzas_elementos


function energias_y_fuerzas_acumuladas_elementos(energ_acumuladas_elementos, inc_energ_elementos, f_acumulado_elementos, inc_u, coords_3D, conectividad_3D, Youngs_elementos)
    N_elementos=length(conectividad_3D[:,1])
    energ_acumuladas_elementos_nuevo=zeros(N_elementos)
    f_acumulado_elementos_nuevo=zeros(24,N_elementos)
    

    for elemento in 1:N_elementos
        u_elemento=zeros(24); inc_u_elemento=zeros(24); coords_vertices_hexaedro=zeros(8,3) # Vector desplazamientos en los nodos del elemento
        for i in 1:8
            nodo=conectividad_3D[elemento,i]
            coords_vertices_hexaedro[i,:]=coords_3D[nodo,1:3]
            inc_u_elemento[3*i-2:3*i]=inc_u[3*nodo-2:3*nodo]
        end

        K_elemento=matriz_K_elemento_hexaedro(coords_vertices_hexaedro, Youngs_elementos[elemento])

        energ_acumuladas_elementos_nuevo[elemento]= energ_acumuladas_elementos[elemento] + inc_energ_elementos[elemento] + f_acumulado_elementos[:,elemento]'*inc_u_elemento
        f_acumulado_elementos_nuevo[:,elemento]= f_acumulado_elementos[:,elemento] + K_elemento*inc_u_elemento

    end

    return energ_acumuladas_elementos_nuevo, f_acumulado_elementos_nuevo

end # function fuerzas_elementos


function resolucion_solido_3D(coords_3D, conectividad_3D, Youngs_elementos, epsilon_voigt_acumulado, gdl_restring_globales, valor_restricciones_globales)

    tipo_invariante="deform_equiv"
   
    ### GENERAR LA MATRIZ DE RIGIDEZ global
    K_global=matriz_K_estructura_3D(coords_3D, conectividad_3D, Youngs_elementos)

    ### OBTENER LOS DESPLAZAMIENTOS EN LOS NODOS
    inc_u,inc_f=solver_desplazamientos_solido(gdl_restring_globales, valor_restricciones_globales, K_global)

    ### INCREMENTOS DE ENERGIA DE LOS ELEMENTOS
    inc_energ_elementos=inc_energias_elementos(inc_u, coords_3D, conectividad_3D, Youngs_elementos)

    ### INVARIANTES DE DEFORMACIONES EN CADA ELEMENTO 
    epsilon_voigt_nuevo,invariante_deform,inc_def_equivalente=invariantes_deformaciones_elementos_hexaedricos(tipo_invariante, coords_3D, conectividad_3D, inc_u, epsilon_voigt_acumulado, Youngs_elementos)
    

    return inc_u, inc_f, inc_energ_elementos, invariante_deform, inc_def_equivalente, epsilon_voigt_nuevo

end # function resolucion_solido_3D 


function resolucion_iterativa_solido(coords_3D, conectividad_3D, Youngs_tangentes_elementos, epsilon_voigt_acumulado, gdl_restring_globales, valor_restricciones_globales, f_acumulado_hexaedros, energias_hexaedros_acumulado)

    Youngs_tangentes_in=copy(Youngs_tangentes_elementos)
    tabla_youngs_tangentes=readdlm(raw"young_tangente_vs_deformacion_solido_homogeneo.txt") #leer el archivo para recoger datos del ensayo
    invariante_deform_max=maximum(tabla_youngs_tangentes[:,2])

    N_elementos=length(Youngs_tangentes_in)
    energias_hexaedros_acumulado_nuevo=[]; fuerzas_hexaedros_acumulado_nuevo=[]
    inc_energia_antiguo=0; inc_energia_antiguo_2=0; n=0
    inc_u=0; inc_f=0; inc_energ_elementos=0; invariante_deform=0; inc_def_equivalente=0; epsilon_voigt_nuevo=0

    iter_max=10; tolerancia1=1e-6; tolerancia2=1e-2
    
    for iter in 1:iter_max
        n=n+1

        ### OBTENER LOS INCREMENTOS DE FUERZA, DESPLAZAMIENTOS...
        inc_u, inc_f, inc_energ_elementos, invariante_deform, inc_def_equivalente, epsilon_voigt_nuevo = resolucion_solido_3D(coords_3D, conectividad_3D, Youngs_tangentes_elementos, epsilon_voigt_acumulado, gdl_restring_globales, valor_restricciones_globales)
        inc_energia=sum(inc_energ_elementos)
        println("iteracion ",iter, ", inc_energia= ", inc_energia, ", deform_maxima= ", maximum(invariante_deform), ", inc_deform_max= ", maximum(inc_def_equivalente)); println()


        ### ENERGIA TOTAL ACUMULADA EN CADA ELEMENTO Y VECTOR DE FUERZAS EN CADA ELEMENTO
        energias_hexaedros_acumulado_nuevo, fuerzas_hexaedros_acumulado_nuevo = energias_y_fuerzas_acumuladas_elementos(energias_hexaedros_acumulado, inc_energ_elementos, f_acumulado_hexaedros, inc_u, coords_3D, conectividad_3D, Youngs_tangentes_elementos)


        ### ACTUALIZAR LOS YOUNGS TANGENTES
        for elemento in 1:N_elementos
            if abs(invariante_deform[elemento])<tabla_youngs_tangentes[1,2]
                Youngs_tangentes_1=tabla_youngs_tangentes[1,1]
            else
                j=1
                while tabla_youngs_tangentes[j,2]-abs(invariante_deform[elemento])<0; 
                    if invariante_deform_max > abs(invariante_deform[elemento])
                        j=j+1; 
                    else
                        println("ERROR: fuera del rango de la curva sigma-epsilon del solido")
                        j=length(tabla_youngs_tangentes[:,2])
                        break
                    end
                end
                Youngs_tangentes_1=(tabla_youngs_tangentes[j,1]+tabla_youngs_tangentes[j-1,1])/2
            end

            Youngs_tangentes_elementos[elemento]=(Youngs_tangentes_1)#+Youngs_tangentes_in[elemento])/2
        end
        # Young_medio=sum(Youngs_tangentes_elementos)/N_elementos
        # Youngs_tangentes_elementos=Young_medio*ones(N_elementos)

        
        if abs(inc_energia-inc_energia_antiguo)/inc_energia < tolerancia1;   break
        else        
            if (abs(inc_energia-inc_energia_antiguo_2)/inc_energia < tolerancia1) && (abs(inc_energia-inc_energia_antiguo)/inc_energia < tolerancia2); break
            else
                inc_energia_antiguo_2=inc_energia_antiguo
                inc_energia_antiguo=inc_energia
            end
        end

    end

    if n==iter_max; println("se ha alcanzado el maximo de iteraciones en el proceso de convergencia"); else; end

    return inc_u, inc_f, inc_energ_elementos, invariante_deform, inc_def_equivalente, epsilon_voigt_nuevo, Youngs_tangentes_elementos, energias_hexaedros_acumulado_nuevo, fuerzas_hexaedros_acumulado_nuevo

end


end # module