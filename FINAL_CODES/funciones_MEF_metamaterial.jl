#MODULO PARA REALIIZAR LOS CÁLCULOS MEF DE DEFORMACIONES, FUERZAS, DESPLAZAMIENTOS, ENERGÍA, etc. DEL METAMATERIAL (ELEMENTOS TIPO BARRA-VIGA)
module funciones_MEF_metamaterial

include(raw"curva_sigma_epsilon_material.jl"); using .curva_sigma_epsilon_material
using LinearAlgebra
using StaticArrays
using Plots
using SparseArrays

export resolucion_iterativa_metamaterial, energia_interna_total_calculo_lineal, condic_contorno, energia_barras_aisladas, generar_diccionarios, tensiones_criticas_pandeo, error_fuerzas, derivadas_nodos_no_fijos, derivadas_nodos_no_fijos_elemento, derivadas_nodos_no_fijos_K


function matriz_K_elemento(coord_inic, coord_fin, Young_elemento, radio)   # En coordenadas globales: elementos viga hermitica
    E = Young_elemento; Pois=Poisson(); G = E/(2*(1+Pois)) #rigidez a cortadura

    A = pi*radio^2   # (area) Por ser vigas circulares
    L = norm(coord_fin[1:3]-coord_inic[1:3])  # longitud del elemento
    #Vol = A*L
    #masa = densidad*Vol
    Iy= 1/4*pi*(radio^2)^2  #momento de inercia respecto al eje y (el eje x es el de la barra)
    Iz= 1/4*pi*(radio^2)^2  #momento de inercia respecto al eje z
    J = Iy+Iz  #modulo a torsion (viga circular)
    

    #Matriz de rigidez de un elemento tipo viga en 3 dimensiones
    K_local =  [E*A/L     0              0           0      0           0        -E*A/L     0           0           0       0           0 ;
                0     12*E*Iz/(L^3)      0           0      0       6*E*Iz/(L^2)    0  -12*E*Iz/(L^3)   0           0       0       6*E*Iz/(L^2); 
                0         0          12*E*Iy/L^3     0  -6*E*Iy/L^2     0           0       0       -12*E*Iy/L^3    0   -6*E*Iy/L^2     0 ;
                0         0              0         G*J/L    0           0           0       0           0         -G*J/L    0           0 ;
                0         0          -6*E*Iy/L^2     0    4*E*Iy/L      0           0       0        6*E*Iy/L^2     0    2*E*Iy/L       0; 
                0      6*E*Iz/L^2        0           0      0        4*E*Iz/L       0   -6*E*Iz/L^2     0           0       0       2*E*Iz/L ;
              -E*A/L      0              0           0      0           0         E*A/L     0           0           0       0           0 ;
                0    -12*E*Iz/L^3        0           0      0       -6*E*Iz/L^2     0    12*E*Iz/L^3    0           0       0       -6*E*Iz/L^2 ;
                0         0          -12*E*Iy/L^3    0    6*E*Iy/L^2    0           0       0       12*E*Iy/L^3     0    6*E*Iy/L^2     0 ;
                0         0              0        -G*J/L    0           0           0       0           0          G*J/L    0           0 ;
                0         0          -6*E*Iy/L^2     0    2*E*Iy/L      0           0       0        6*E*Iy/L^2     0     4*E*Iy/L      0; 
                0    6*E*Iz/L^2          0           0      0         2*E*Iz/L      0    -6*E*Iz/L^2    0           0       0        4*E*Iz/L]
    
    # Pasar de coordenadas locales a globales (matriz de cambio de base)
    R_12_12= ejes_locales_a_globales(coord_inic,coord_fin) 

    # matriz de rigidez del elemento en globales
    K_elemento =R_12_12*(K_local*R_12_12')  

    return  K_elemento

end #function matriz_K_elemento


function ejes_locales_a_globales(coord_inic,coord_fin) # Matriz cambio de base para pasar de coordenadas locales a globales
    L = norm(coord_fin[1:3]-coord_inic[1:3])  # longitud del elemento

    # eje x unitario de las cooordenadas locales y ejes Y, Z (elegidos arbitrariamiente pero perpendiculares)
    eje_x_local= (coord_fin[1:3]-coord_inic[1:3])/L   
    if eje_x_local[1]!=0 || eje_x_local[2]!=0  # eje prependicular al X local (y unitario) --> hay que definir de dos formas al vector Y local para evitar que sea nulo en algunos caso
        eje_y_local= [eje_x_local[2],-eje_x_local[1],0];  eje_y_local=eje_y_local/norm(eje_y_local)  
    else
        eje_y_local=[0,-eje_x_local[3],eje_x_local[2]];  eje_y_local=eje_y_local/norm(eje_y_local)
    end
    eje_z_local= cross(eje_x_local,eje_y_local)

    eje_x_global=[1,0,0]; eje_y_global=[0,1,0]; eje_z_global=[0,0,1];

    #matriz de cambio de coordenadas de 3x3 (local a global) y de 12x12 (para el elemento completo con giros y desplazamientos)
    R_3_3 = [dot(eje_x_local,eje_x_global) dot(eje_y_local,eje_x_global) dot(eje_z_local,eje_x_global);
             dot(eje_x_local,eje_y_global) dot(eje_y_local,eje_y_global) dot(eje_z_local,eje_y_global);
             dot(eje_x_local,eje_z_global) dot(eje_y_local,eje_z_global) dot(eje_z_local,eje_z_global)]

    R_12_12 = zeros(12,12)
    R_12_12[1:3,1:3] = R_3_3
    R_12_12[4:6,4:6] = R_3_3
    R_12_12[7:9,7:9] = R_3_3
    R_12_12[10:12,10:12] = R_3_3

    return R_12_12

end #function desplaz_globales_a_locales


function matriz_K_estructura(conectividad_zona, coords_zona_metam, Youngs_tangentes_barras, radio)
    N_elementos=length(conectividad_zona[:,1])
    #obtenermos los nodos a partir de la matriz de conectividad_zona    
    N_gdl_zona=6*length(coords_zona_metam[:,1])
    K_global=zeros(N_gdl_zona,N_gdl_zona)


    for n in 1:N_elementos
        
        nodo_inic_local=findall(in(conectividad_zona[n,1]),coords_zona_metam[:,4])[1] # la numeracion local da valores del 1,...,N_nodos a cada nodo 
        nodo_fin_local=findall(in(conectividad_zona[n,2]),coords_zona_metam[:,4])[1]
        coords_nodo_inic=coords_zona_metam[nodo_inic_local,1:3]    #coordenadas del nodo inicial
        coords_nodo_fin=coords_zona_metam[nodo_fin_local,1:3]    #coordenadas del nodo inicial

        gdl_elemento=vcat(collect(6*nodo_inic_local-5:6*nodo_inic_local),collect(6*nodo_fin_local-5:6*nodo_fin_local))

        elemento_metam=conectividad_zona[n,3]
        #println(Youngs_tangentes_barras)
        Young_elemento=Youngs_tangentes_barras[elemento_metam]

        #matriz K del elemento asiciado a los nodos anteriores
        K_elemento = matriz_K_elemento(coords_nodo_inic, coords_nodo_fin, Young_elemento, radio)
     
        #ensamblaje de la matriz del elemento en la matriz completa
        K_global[gdl_elemento,gdl_elemento]=K_global[gdl_elemento,gdl_elemento] + K_elemento   

    end  

    return K_global

end #function matriz_K_estructura


function condic_contorno(punto_restriccion, gdl_restriccion, punto_desplazamiento, valor_desplazamiento, coords) # condic de contorno en desplazamientos
    ## DETERMINA LOS GDL CON DESPLAZAMIENTO IMPUESTO Y CON RESTRICCIONES
    N_restricciones= length(punto_restriccion[:,1])
    N_desplazamientos= length(punto_desplazamiento[:,1])
    N_nodos= length(coords[:,1])
    gdl_desplazamiento=valor_desplazamiento./valor_desplazamiento;  gdl_desplazamiento=replace(gdl_desplazamiento, NaN=>0)
    gdl_restring_globales=[];  gdl_desplaz_globales=[];  valor_desplaz_globales=[]#zeros(6*N_desplazamientos,2); 
    
    epsilon=1e-8
    coordenadas_x=sort(unique(coords[:,1])); coordenadas_y=sort(unique(coords[:,2])); coordenadas_z=sort(unique(coords[:,3]))
    dimension_x_celdilla=coordenadas_x[2]-coordenadas_x[1];
    dimension_y_celdilla=coordenadas_y[2]-coordenadas_y[1];
    dimension_z_celdilla=coordenadas_z[2]-coordenadas_z[1];


    # RESTRICCIONES QUE IMPIDEN EL MOVIMIENTO
    # se recorren todos los nodos y se comprueba se se situan en una coordenada donde haya impuesta una condicion de contorno
    for n in 1:N_restricciones

        a=length(gdl_restring_globales) # esto se hace para ver si la condicion de contorno no coincide exactamente con ningun nodo
        for m in 1:N_nodos 

            gdlx=Int64(gdl_restriccion[n,1]*(6*m-5))
            gdly=Int64(gdl_restriccion[n,2]*(6*m-4))
            gdlz=Int64(gdl_restriccion[n,3]*(6*m-3))
            gdlthetax=Int64(gdl_restriccion[n,4]*(6*m-2))
            gdlthetay=Int64(gdl_restriccion[n,5]*(6*m-1))
            gdlthetaz=Int64(gdl_restriccion[n,6]*(6*m))
            gdl_restring_globales_nodo=[gdlx,gdly,gdlz,gdlthetax,gdlthetay,gdlthetaz]

            if punto_restriccion[n,1]=="all" || abs(coords[m,1]-punto_restriccion[n,1]) < epsilon
                if punto_restriccion[n,2]=="all"  || abs(coords[m,2]-punto_restriccion[n,2]) < epsilon
                    if punto_restriccion[n,3]=="all" || abs(coords[m,3]-punto_restriccion[n,3]) < epsilon
                        append!(gdl_restring_globales,gdl_restring_globales_nodo)
                        #println(coords[m,:])
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

                gdlx=Int64(gdl_restriccion[n,1]*(6*m-5))
                gdly=Int64(gdl_restriccion[n,2]*(6*m-4))
                gdlz=Int64(gdl_restriccion[n,3]*(6*m-3))
                gdlthetax=Int64(gdl_restriccion[n,4]*(6*m-2))
                gdlthetay=Int64(gdl_restriccion[n,5]*(6*m-1))
                gdlthetaz=Int64(gdl_restriccion[n,6]*(6*m))
                gdl_restring_globales_nodo=[gdlx,gdly,gdlz,gdlthetax,gdlthetay,gdlthetaz]

                if punto_restriccion[n,1]=="all" || abs(coords[m,1]-punto_restriccion[n,1]) < dimension_x_celdilla
                    if punto_restriccion[n,2]=="all"  || abs(coords[m,2]-punto_restriccion[n,2]) < dimension_y_celdilla
                        if punto_restriccion[n,3]=="all" || abs(coords[m,3]-punto_restriccion[n,3]) < dimension_z_celdilla
                            append!(gdl_restring_globales,gdl_restring_globales_nodo)
                            #println(coords[m,:])
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

            gdlx=Int64(gdl_desplazamiento[n,1]*(6*m-5));           valor_desplaz_x=valor_desplazamiento[n,1]
            gdly=Int64(gdl_desplazamiento[n,2]*(6*m-4));           valor_desplaz_y=valor_desplazamiento[n,2]
            gdlz=Int64(gdl_desplazamiento[n,3]*(6*m-3));           valor_desplaz_z=valor_desplazamiento[n,3]
            gdlthetax=Int64(gdl_desplazamiento[n,4]*(6*m-2));      valor_theta_x=valor_desplazamiento[n,4]
            gdlthetay=Int64(gdl_desplazamiento[n,5]*(6*m-1));      valor_theta_y=valor_desplazamiento[n,5]
            gdlthetaz=Int64(gdl_desplazamiento[n,6]*(6*m));        valor_theta_z=valor_desplazamiento[n,6]
            gdl_desplaz_nodo=[gdlx,gdly,gdlz,gdlthetax,gdlthetay,gdlthetaz]
            valor_desplaz_nodo=[valor_desplaz_x,valor_desplaz_y,valor_desplaz_z,valor_theta_x,valor_theta_y,valor_theta_z]

            if punto_desplazamiento[n,1]=="all" || abs(coords[m,1]-punto_desplazamiento[n,1]) < epsilon
                if punto_desplazamiento[n,2]=="all"  || abs(coords[m,2]-punto_desplazamiento[n,2]) < epsilon
                    if punto_desplazamiento[n,3]=="all" || abs(coords[m,3]-punto_desplazamiento[n,3]) < epsilon
                        append!(gdl_desplaz_globales, gdl_desplaz_nodo)
                        append!(valor_desplaz_globales, valor_desplaz_nodo)
                        #i=i+1
                        #println(coords[m,:])
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
            # Esta aproximacion al nodo mas cercano será mas cierta cuanto mas celdillas haya
            for m in 1:N_nodos 

                gdlx=Int64(gdl_desplazamiento[n,1]*(6*m-5));           valor_desplaz_x=valor_desplazamiento[n,1]
                gdly=Int64(gdl_desplazamiento[n,2]*(6*m-4));           valor_desplaz_y=valor_desplazamiento[n,2]
                gdlz=Int64(gdl_desplazamiento[n,3]*(6*m-3));           valor_desplaz_z=valor_desplazamiento[n,3]
                gdlthetax=Int64(gdl_desplazamiento[n,4]*(6*m-2));      valor_theta_x=valor_desplazamiento[n,4]
                gdlthetay=Int64(gdl_desplazamiento[n,5]*(6*m-1));      valor_theta_y=valor_desplazamiento[n,5]
                gdlthetaz=Int64(gdl_desplazamiento[n,6]*(6*m));        valor_theta_z=valor_desplazamiento[n,6]
                gdl_desplaz_nodo=[gdlx,gdly,gdlz,gdlthetax,gdlthetay,gdlthetaz]
                valor_desplaz_nodo=[valor_desplaz_x,valor_desplaz_y,valor_desplaz_z,valor_theta_x,valor_theta_y,valor_theta_z]
    
                if punto_desplazamiento[n,1]=="all" || abs(coords[m,1]-punto_desplazamiento[n,1]) < dimension_x_celdilla
                    if punto_desplazamiento[n,2]=="all"  || abs(coords[m,2]-punto_desplazamiento[n,2]) < dimension_y_celdilla
                        if punto_desplazamiento[n,3]=="all" || abs(coords[m,3]-punto_desplazamiento[n,3]) < dimension_z_celdilla
                            append!(gdl_desplaz_globales, gdl_desplaz_nodo)
                            append!(valor_desplaz_globales, valor_desplaz_nodo)
                            #i=i+1
                            #println(coords[m,:])
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
    #gdl_desplaz_globales=hcat(gdl_desplaz_globales,valor_desplaz_globales)
    
    # Hay que eliminar los gdl de los nodos a los que se les ha impuesto una condicion porque no en todos los gdl tienen una
    # lo hacamos para las restricciones que impiden el movimiento
    gdl_no_nulo=findall(!in(0),gdl_restring_globales); gdl_restring_globales=gdl_restring_globales[gdl_no_nulo]
    # lo hacemos para las restricciones con un desplazamiento asociado
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


    #ver si algun nodo tiene impuestas restricciones y desplazamientos a la vez
    if findall(in(gdl_desplaz_globales[:,1]),gdl_restring_globales)!=[]
        println("un mismo nodo tiene impuestos una restriccion y un desplazamiento a la vez")
    else
    end   


    # OBTENER TODOS LOS GDL CON ALGUNA RESTRICCION
    N_restricciones_movimiento_impedido=length(gdl_restring_globales)
    append!(gdl_restring_globales, gdl_desplaz_globales) #todos los nodos con ALGUNA condicion de contorno (y ordenados)
    valor_restricciones_globales=append!(zeros(N_restricciones_movimiento_impedido), valor_desplaz_globales)


    # Ver si a algun grado de libertad se repite porque se le han aplicado 2 condiciones de contorno diferentes
    if length(gdl_restring_globales) > length(unique(gdl_restring_globales))
        println("2 condiciones distintas para un mismo gdl")
    else
    end

    return gdl_restring_globales, valor_restricciones_globales

end #function condic_contorno


function solver_desplazamientos(gdl_restring, valor_restriccion, K_global) #obtener el vector u de desplazamientos globales
    
    gdl=length(K_global[:,1]); N_nodos=Int64(gdl/6)   #numero de grados de libertad y de nodos
    gdl_totales=collect(1:gdl)  #vector con los grados de libertad
    gdl_libres=findall(!in(gdl_restring),gdl_totales) #gdl sin restricciones   
    u=zeros(gdl); f=zeros(gdl) #desplazamiento y fuerza

    # Resolucion por condensacion estatica
    u[gdl_restring]=valor_restriccion #asigno la condicion de contorno donde corresponda
    #$$$$
    dimension_problema=length(gdl_libres);println("Dimension del problema de metamaterial: ", dimension_problema);
    dimension_limite=1000
    iter_max=10000

    if dimension_problema <= dimension_limite # con dimensiones pequeñas se resuelve la inversa numericamente

        u[gdl_libres]=K_global[gdl_libres,gdl_libres]\(-K_global[gdl_libres,gdl_restring]*u[gdl_restring]) #desplazamiento de los puntos sin ninguna restriccion
    
    else # con dimensiones grandes se resuelve la inversa iterativamente

        A=sparse(K_global[gdl_libres,gdl_libres]); b=-K_global[gdl_libres,gdl_restring]*u[gdl_restring]

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

    f=K_global*u #f[gdl_restring]=K_global[gdl_restring,gdl_libres]*u[gdl_libres]+K_global[gdl_restring,gdl_restring]*u[gdl_restring]
    #$$$$
    #println("Problema de metamaterial resuelto")

    # separar desplazamientos y giros asi como fuerzas y momentos
    # giros=zeros(3*N_nodos); desplazamientos=zeros(3*N_nodos); fuerzas=zeros(3*N_nodos); momentos=zeros(3*N_nodos)
    # n=1
    # for i in 1:6:gdl
    #     desplazamientos[n:n+2]=u[i:i+2]
    #     giros[n:n+2]=u[i+3:i+5]
    #     fuerzas[n:n+2]=f[i:i+2]
    #     momentos[n:n+2]=f[i+3:i+5]
    #     n=n+3
    # end

    return u,f

end #function solver_desplazamientos


function deformaciones_ptos_int(inc_u_dict, eps_x_ptos_int_acumulado, eps_xy_ptos_int_acumulado, eps_xz_ptos_int_acumulado, coords_zona_metam, conectividad_zona, radio) # deformacion en los puntos de integracion en todos los elementos (elemento viga hermitica)
    #el punto elegido dependerá de si queremos usar la deformacion maxima, media o en los puntos de ...
    #...integracion para representar la del elemento completo (cogiendo la máxima se está del lado de la seguridad)

    N_elementos=length(conectividad_zona[:,1])
    chi_ptos_int=[-1/sqrt(3),1/sqrt(3)]   #coordenada naturales de los puntos de integración
    eps_x_ptos_int=Dict(); eps_xy_ptos_int=Dict(); eps_xz_ptos_int=Dict(); # deformaciones en los puntos de integracion

    #$$$$
    a=[];b=[];c=[]

    #incrementos de desplazamientos en ejes locales para cada elemento y deformacion
    for n in 1:N_elementos
        u_elemento_glob=zeros(12) #vector de desplazamiento del elemento en globales
        elemento_metam=conectividad_zona[n,3]

        # obtener los nodos que forman el elemento barra
        nodo_inic=conectividad_zona[n,1] #nodo inicial de un elemento
        indice_nodo_inic=findall(in(nodo_inic),coords_zona_metam[:,4])[1]    #coordenadas del nodo inicial
        coords_nodo_inic=coords_zona_metam[indice_nodo_inic,1:3]
        nodo_fin=conectividad_zona[n,2] #nodo final de un elemento
        indice_nodo_fin=findall(in(nodo_fin),coords_zona_metam[:,4])[1]    #coordenadas del nodo inicial
        coords_nodo_fin=coords_zona_metam[indice_nodo_fin,1:3]

        
        # OBTENER LOS INCREMENTOS DE DESPLAZAMIENTOS DEL ELEMENTO
        u_elemento_glob[1:6]=inc_u_dict[nodo_inic] #[gdl_inic_1_traducido:gdl_inic_1_traducido+5]
        u_elemento_glob[7:12]=inc_u_dict[nodo_fin] #[gdl_fin_1_traducido:gdl_fin_1_traducido+5]


        R_12_12=ejes_locales_a_globales(coords_nodo_inic,coords_nodo_fin) #matriz para pasar de coordenadas locales a globales

        u_elemento_loc=R_12_12'*u_elemento_glob
        u_x1=u_elemento_loc[1];             u_x2=u_elemento_loc[7];
        u_y1=u_elemento_loc[2];             u_y2=u_elemento_loc[8];
        u_z1=u_elemento_loc[3];             u_z2=u_elemento_loc[9];
        theta_x1=u_elemento_loc[4];         theta_x2=u_elemento_loc[10];
        theta_y1=u_elemento_loc[5];         theta_y2=u_elemento_loc[11];
        theta_z1=u_elemento_loc[6];         theta_z2=u_elemento_loc[12];
        #println(elemento);println("u_glob= ",u_elemento_glob);println("u_loc= ",u_elemento_loc)


        # Jacobiano del elemento
        L = norm(coords_nodo_fin-coords_nodo_inic)  # longitud del elemento
        Jac=L/2; inv_Jac=1/Jac

        y=radio; z=radio; r=radio/2             #para obtener el valor máximo de deformacion en los puntos de integracion

        #$$$$
        inc_eps_x_axil=0;

        eps_x_ptos_int_acumulad=eps_x_ptos_int_acumulado[elemento_metam]; eps_xy_ptos_int_acumulad=eps_xy_ptos_int_acumulado[elemento_metam]; eps_xz_ptos_int_acumulad=eps_xz_ptos_int_acumulado[elemento_metam]
        eps_x_ptos_int[elemento_metam]=[]; eps_xy_ptos_int[elemento_metam]=[]; eps_xz_ptos_int[elemento_metam]=[]

        for i in 1:length(chi_ptos_int)

            ## INCREMENTO DEFORMACION epsilon_x POR PARTE AXIL (2 puntos de integracion)--> salvo que sea nulo, el axil suele ser el dominante
            N=1/2*[1-chi_ptos_int[i],1+chi_ptos_int[i]]
            B=inv_Jac*1/2*[-1,1]
            inc_eps_x_axil=dot(B,[u_x1,u_x2])

            ### INCREMENTO DEFORMACION epsilon_x POR PARTE DE FLEXION EN Y y en Z LOCAL a una distancia del eje (radio) (viga hermitica)
            N_herm_1=-0.75*chi_ptos_int[i]+0.25*chi_ptos_int[i]^3+0.5
            N_herm_2=-0.25*chi_ptos_int[i]-0.25*chi_ptos_int[i]^2+0.25*chi_ptos_int[i]^3+0.25
            N_herm_3=0.75*chi_ptos_int[i]-0.25*chi_ptos_int[i]^3+0.5
            N_herm_4=-0.25*chi_ptos_int[i]+0.25*chi_ptos_int[i]^2+0.25*chi_ptos_int[i]^3-0.25
            N_herm=[N_herm_1,Jac*N_herm_2,N_herm_3,Jac*N_herm_4]
            B_herm=inv_Jac^2*[6/4*chi_ptos_int[i], Jac*(-2+6*chi_ptos_int[i])/4, -6/4*chi_ptos_int[i], Jac*(2+6*chi_ptos_int[i])/4]

            # Hay que determinar el punto de la seccion de los ptos de integracion (coordenadas Y y Z) donde la deformacion es maxima (será el mas critico)
            inc_eps_x_flex_y_aux=dot(B_herm,[u_y1,theta_z1,u_y2,theta_z2])
            inc_eps_x_flex_z_aux=dot(B_herm,[u_z1,theta_y1,u_z2,theta_y2])

            
            alpha=atan(-inc_eps_x_flex_z_aux,inc_eps_x_flex_y_aux) #resolviendo: maximo(-radio*cos(alpha)*inc_eps_x_flex_y_aux+radio*sin(alpha)*inc_eps_x_flex_z_aux)

            inc_eps_x_flex_y=r*cos(alpha)*inc_eps_x_flex_y_aux
            inc_eps_x_flex_z=r*sin(alpha)*inc_eps_x_flex_z_aux
            inc_eps_x_flex=-inc_eps_x_flex_y + inc_eps_x_flex_z


            ### Se debe tener en cuenta si eps_axil acumulado en el elemento es positivo o negativo a la hora de seleccionar el inc_epsilon de flexion que hace mas critica a la viga
            if inc_eps_x_axil > 0
                push!(eps_x_ptos_int[elemento_metam], eps_x_ptos_int_acumulad[i]+inc_eps_x_axil+abs(inc_eps_x_flex))
            elseif inc_eps_x_axil < 0
                push!(eps_x_ptos_int[elemento_metam], eps_x_ptos_int_acumulad[i]+inc_eps_x_axil-abs(inc_eps_x_flex))
            else #si el axil es cero se coge el positivo (de traccion) (criterio arbitrario)
                push!(eps_x_ptos_int[elemento_metam], eps_x_ptos_int_acumulad[i]+abs(inc_eps_x_flex))
            end


            ### INCREMENTOS DE DEFORMACION gamma_x_r a una distancia del eje (r) POR EL GIRO ENTORNO AL EJE X
            inc_gamma_x_r=r*dot(B,[theta_x1,theta_x2])

            ### a partir del gamma_x_r, obtenemos gamma_xy y gamma_xz en el punto donde es máximo el inc_epsilon_x --> eso se determina por el alpha
            inc_gamma_xy=-inc_gamma_x_r*sin(alpha) #el signo no importa porque se va a usar para los invariantes y alli se usa al cuadrado
            inc_gamma_xz=inc_gamma_x_r*cos(alpha)
            push!(eps_xy_ptos_int[elemento_metam], eps_xy_ptos_int_acumulad[i]+1/2*inc_gamma_xy)
            push!(eps_xz_ptos_int[elemento_metam], eps_xz_ptos_int_acumulad[i]+1/2*inc_gamma_xz)
            
            #println(elemento_metam," ",inc_eps_x_axil," ",inc_eps_x_flex)

        end
        
        #$$$$$$$
        #push!(a,eps_x_axil);push!(b,eps_x_ptos_int[elemento_metam][1]);push!(c,eps_x_ptos_int[elemento_metam][2]) 
    end
    #$$$$$$$
    #display(plot(a,label="eps_axil"));display(plot!(b,label="eps_x_1"));display(plot!(c,label="eps_x_2"))

    return eps_x_ptos_int, eps_xy_ptos_int, eps_xz_ptos_int

end # function deformaciones_ptos_int


function invariantes_deformaciones(tipo_invariante, eps_x_ptos_int, eps_xy_ptos_int, eps_xz_ptos_int, conectividad_zona) # segun el criterio elegido para comparar con el caso unidimensional y los puntos que se quieran coger para comparar, esta funcion cambiará
    ### Criterios relacionados con DEFORMACION MÁXIMA: debe estar en uno de los dos extremos por ser la interpolacion lineal
    # obtenemos los 3 invariantes maximos en cada elemento
    N_elementos=length(conectividad_zona[:,1])
    invariantes_elementos=Dict() #; invariante_2_elemento=Dict(); invariante_3_elemento=Dict(); deformacion_equivalente_elemento=Dict()
    plots_elementos=Dict()
    poisson=Poisson()

    chi_min=-sqrt(3); chi_max=sqrt(3); # Nueva base con los puntos de integración como nodos (homotecia)
    N_puntos_invariantes_barra=5 # numero de puntos dentro de la barra donde estudiamos los invariantes para hacer la media
    
    for n in 1:N_elementos
        elemento_metam=conectividad_zona[n,3]
        invariantes_barra=[]#; invariantes_2=[]; invariantes_3=[]; deformaciones_equivalentes=[]  # invariantes en distintos puntos de un elemento

        for chi in range(chi_min,chi_max,length=N_puntos_invariantes_barra)  # se calculan los invariantes en varios puntos para luego coger el que dé mayores valores; mas puntos = mas precision para coger el maximo 
            N=1/2*[1-chi,1+chi]
            eps_x=dot(N,eps_x_ptos_int[elemento_metam])
            eps_xy=dot(N,eps_xy_ptos_int[elemento_metam])
            eps_xz=dot(N,eps_xz_ptos_int[elemento_metam])
            eps_y=-poisson*eps_x; eps_z=-poisson*eps_x #estimamos el eps_y y el eps_x

            signo=sign(eps_x)

            ### INVARIANTES ESCALADOS (PARA QUE SEAN DEL MISMO ORDEN DE MAGNITUD)
            # este "if" se hace para no hevitar calculos innecesarios
            if tipo_invariante=="inv_1"
                invariante_1=eps_x + eps_y + eps_z
                append!(invariantes_barra,invariante_1)

            elseif tipo_invariante=="inv_2"
                invariante_2=eps_x*eps_y + eps_x*eps_z + eps_z*eps_y - eps_xy*eps_xy - eps_xz*eps_xz
                append!(invariantes_barra,signo*sqrt(abs(invariante_2)))

            elseif tipo_invariante=="inv_3"
                tensor_epsilon=SMatrix{3, 3, Float64,9}([eps_x  eps_xy eps_xz;
                                                         eps_xy  eps_y   0;
                                                         eps_xz   0    eps_z])
                invariante_3=det(tensor_epsilon)
                append!(invariantes_barra,signo*cbrt(abs(invariante_3)))

            elseif tipo_invariante=="deform_equiv"
                deformacion_equivalente=sqrt(1/((1+poisson)*(1-2*poisson))*((1-poisson)*(eps_x^2+eps_y^2+eps_z^2)+2*poisson*(eps_x*eps_y + eps_x*eps_z + eps_z*eps_y))+4/(2*(1+poisson))*(eps_xy^2+eps_xz^2))
                append!(invariantes_barra,signo*deformacion_equivalente)
            
            elseif tipo_invariante=="deform_equiv_2"
                deformacion_equivalente_2=1/(sqrt(2)*(1+poisson))*sqrt((eps_x-eps_y)^2+(eps_y-eps_z)^2+(eps_z-eps_x)^2+6*(eps_xy^2+eps_xz^2))
                append!(invariantes_barra,signo*deformacion_equivalente_2)

            else
                println("No se ha introducido bien el tipo de invariante")
            end

            println(eps_x)    
        end
        println(elemento_metam, " ", invariantes_barra)

        # Obtener los 3 invariantes máximos y minimos en un elemento y así para todos los elementos (para tener en cuenta la traccion y compresion)
        invariante_max=maximum(invariantes_barra);  invariante_min=minimum(invariantes_barra)
        invariante_medio=sum(abs.(invariantes_barra))/N_puntos_invariantes_barra
        # invariante_2_max=maximum(invariantes_2);  invariante_2_min=minimum(invariantes_2)
        # invariante_3_max=maximum(invariantes_3);  invariante_3_min=minimum(invariantes_3) 
        # deform_equiv_max=maximum(deformaciones_equivalentes);  deform_equiv_min=minimum(deformaciones_equivalentes);
        
        #$$$$
        #para ver la distribucion dentro de invariantes elemento
        #plots_elementos[elemento_metam]=plot(invariantes_1,title="invariantes",label="inv_1");plot!(invariantes_2,label="inv_2");plot!(invariantes_3,label="inv_3");plot!(deformaciones_equivalentes,label="def_equiv");xlabel!("elemento");ylabel!("invariantes")
        #readline()

        # OBTENER EL MAXIMO VALOR ABSOLUTO DE LOS INVARIANTES (CON SU SIGNO) (Asi con los 3 invariantes)
        if abs(invariante_max) >= abs(invariante_min);   invariantes_elementos[elemento_metam]=invariante_medio
        else;                                            invariantes_elementos[elemento_metam]=-invariante_medio
        end

        # if abs(invariante_2_max) >= abs(invariante_2_min); invariante_2_elemento[elemento_metam]=invariante_2_max
        # else;                                              invariante_2_elemento[elemento_metam]=invariante_2_min
        # end

        # if abs(invariante_3_max) >= abs(invariante_3_min); invariante_3_elemento[elemento_metam]=invariante_3_max
        # else;                                              invariante_3_elemento[elemento_metam]=invariante_3_min
        # end

        # if abs(deform_equiv_max) >= abs(deform_equiv_min); deformacion_equivalente_elemento[elemento_metam]=deform_equiv_max
        # else;                                              deformacion_equivalente_elemento[elemento_metam]=deform_equiv_min
        # end
    end
    
    #$$$$
    ### PLOTEAR LOS INVARIANTES PARA COMPARARLOS EN LOS DISTINTOS ELEMENTOS
    # a=[]; b=[];c=[];d=[]
    # for i in conectividad_zona[:,3]
    #     push!(a,invariante_1_elemento[i]); push!(b,invariante_2_elemento[i]); push!(c,invariante_3_elemento[i]); push!(d,deformacion_equivalente_elemento[i]);
    #     println(invariante_1_elemento[i]," ",invariante_2_elemento[i]," ",invariante_3_elemento[i]," ",deformacion_equivalente_elemento[i]);
    # end
    # display(plot!(a,title="invariantes",label="inv_1"));display(plot!(b,label="inv_2"));display(plot!(c,label="inv_3"));display(plot!(d,label="def_equiv"));display(xlabel!("elementos"));display(ylabel!("invariantes"))


    return  invariantes_elementos #, invariante_2_elemento, invariante_3_elemento, deformacion_equivalente_elemento, plots_elementos

end  #function invariantes_deformaciones


# La energia interna del conjunto se puede calcular solo como 1/2*u_global'*K_global*u_global
function energia_interna_total_calculo_lineal(conectividad_zona, coords_zona_metam, u_reducido_dict, Youngs_in, radio) #calcula la energia interna del conjunto de todas las barras como suma de la de cada barra
    N_elementos=length(conectividad_zona[:,1])
    energ_interna_elemento=Dict()


    #desplazamientos en ejes locales para cada elemento y deformacion
    for n in 1:N_elementos

        u_elemento_glob=zeros(12) #vector de desplazamiento del elemento en globales
        elemento_metam=conectividad_zona[n,3]

        # obtener los nodos que forman el elemento barra
        nodo_inic=conectividad_zona[n,1] #nodo inicial de un elemento
        indice_nodo_inic=findall(in(nodo_inic),coords_zona_metam[:,4])    #coordenadas del nodo inicial
        coords_nodo_inic=coords_zona_metam[indice_nodo_inic,1:3]
        nodo_fin=conectividad_zona[n,2] #nodo final de un elemento
        indice_nodo_fin=findall(in(nodo_fin),coords_zona_metam[:,4])    #coordenadas del nodo inicial
        coords_nodo_fin=coords_zona_metam[indice_nodo_fin,1:3]
    
        u_elemento_glob[1:6]=u_reducido_dict[nodo_inic] #[gdl_inic_1_traducido:gdl_inic_1_traducido+5]
        u_elemento_glob[7:12]=u_reducido_dict[nodo_fin] #[gdl_fin_1_traducido:gdl_fin_1_traducido+5]


        ### modulo de Young de cada barra
        Young_elemento=Youngs_in[elemento_metam]

        K_elemento=matriz_K_elemento(coords_nodo_inic, coords_nodo_fin, Young_elemento, radio)

        energ_interna_elemento[elemento_metam]=1/2*u_elemento_glob'*(K_elemento*u_elemento_glob)  # Suponiendo que se puede calcular como en un elemento lineal: U=1/2*u*K*u
    
    end

    energia_interna_total=sum(values(energ_interna_elemento))

    return energ_interna_elemento

end #function energia_interna_total


function energia_interna_total_curva_tension_deformacion(energias_unitarias_elementos, conectividad_zona, coords_zona_metam, radio) # obtener la energia interna de la zona a partir de la integral de la curva tension-deformacion
    N_elementos=length(conectividad_zona[:,1])
    energias_elementos=Dict()

    for n in 1:N_elementos
        elemento=conectividad_zona[n,3]

        nodo_inic=conectividad_zona[n,1] #nodo inicial de un elemento
        indice_nodo_inic=findall(in(nodo_inic),coords_zona_metam[:,4])    #coordenadas del nodo inicial
        coords_nodo_inic=coords_zona_metam[indice_nodo_inic,1:3]
        nodo_fin=conectividad_zona[n,2] #nodo final de un elemento
        indice_nodo_fin=findall(in(nodo_fin),coords_zona_metam[:,4])    #coordenadas del nodo inicial
        coords_nodo_fin=coords_zona_metam[indice_nodo_fin,1:3]

        volumen_elemento=norm(coords_nodo_fin-coords_nodo_inic)*pi*radio^2

        # la energia de un elemento será la energia por unidad de volumen multiplicada por el volumen
        energias_elementos[elemento] = energias_unitarias_elementos[elemento]*volumen_elemento

    end

    return energias_elementos

end #function energia_interna_total_no_lineal



function resolucion_iterativa_metamaterial(coords_zona_metam, conectividad_zona, gdl_restring_globales_metam, valor_restricciones_globales_metam, Youngs_tangentes_in, tensiones_pandeo ,eps_x_ptos_int_acumulado, eps_xy_ptos_int_acumulado, eps_xz_ptos_int_acumulado, energias_acumuladas_elementos, f_acumulada_elementos, invariante_def, radio) # Proceso iterativo de resolucion del problema (necesario cuando se supera el limite elastico)
    
    Youngs_tangentes_barras_actual =copy(Youngs_tangentes_in);    Youngs_tangentes_barras =copy(Youngs_tangentes_in) 
    tipo_invariante = "deform_equiv"  # el invariante escogido para hacer cálculos 

    iter_max=10; tolerancia1=1e-6; tolerancia2=1e-2

    inc_u_dict=0; coords_nodos_elemento=Dict()
    f_acumulada_elementos_nuevo=Dict(); energias_acumuladas_elementos_nuevo=Dict()
    inc_energia=0; inc_energia_antiguo=0; inc_energia_antiguo_2=0; energia_plot=0
    inc_u=0; inc_f=0; inc_energias_elementos=Dict(); eps_x_ptos_int_nuevo=0; eps_xy_ptos_int_nuevo=0; eps_xz_ptos_int_nuevo=0; invariante_def_nuevo=Dict()
    
    n=0

    for iter in 1:iter_max; n=n+1

        ### PARA OBTENER LA PENDIENTE EN EL METODO PREDICTOR-CORRECTOR
        for elemento in conectividad_zona[:,3]
          Youngs_tangentes_barras[elemento]=(Youngs_tangentes_barras_actual[elemento]+Youngs_tangentes_in[elemento])/2
        end          


        ### NUEVA MATRIZ DE RIGIDEZ DE LA ESTRUCTURA DE METAMATERIAL
        K_zona =matriz_K_estructura(conectividad_zona, coords_zona_metam, Youngs_tangentes_barras, radio)


        ### OBTENER LOS INCREMENTOS DE DESPLAZAMIENTOS, DE FUERZAS,...:
        inc_u,inc_f=solver_desplazamientos(gdl_restring_globales_metam, valor_restricciones_globales_metam, K_zona)
        inc_u_dict=generar_diccionario_u(inc_u, coords_zona_metam)


        ### OBTENER EL INCREMENTO DE ENERGIA DE CADA UNA DE LAS BARRAS Y SU ENERGIA ACUMULADA
        for i in 1:length(conectividad_zona[:,3])
            nodo_1, nodo_2, elemento=conectividad_zona[i,:]
            if haskey(coords_nodos_elemento, elemento)==false
                indice_1=findall(in(nodo_1),coords_zona_metam[:,4]); coords_nodo_1=coords_zona_metam[indice_1,1:3]
                indice_2=findall(in(nodo_2),coords_zona_metam[:,4]); coords_nodo_2=coords_zona_metam[indice_2,1:3]
                
                coords_nodos_elemento[elemento]=[coords_nodo_1;  coords_nodo_2]
            else
                coords_nodo_1=coords_nodos_elemento[elemento][1,1:3]
                coords_nodo_2=coords_nodos_elemento[elemento][2,1:3]
            end

            inc_u_elemento=vcat(inc_u_dict[nodo_1],inc_u_dict[nodo_2])
            K_elemento = matriz_K_elemento(coords_nodo_1, coords_nodo_2, Youngs_tangentes_barras[elemento], radio)

            inc_f_elemento=K_elemento*inc_u_elemento
            inc_energia_elemento=1/2*inc_u_elemento'*inc_f_elemento
            inc_energias_elementos[elemento]=inc_energia_elemento

            L=norm(coords_nodo_1-coords_nodo_2)
            invariante_def_nuevo[elemento]=invariante_def[elemento]+sqrt(inc_energia_elemento/(1/2*Youngs_tangentes_barras[elemento]*pi*radio^2*L))
            #println(elemento," ",invariante_def[elemento]," ",invariante_def_nuevo[elemento], " ", inc_energias_elementos[elemento])
            #f_acumulada_elementos_nuevo[elemento]=f_acumulada_elementos[elemento] + inc_f_elemento
            #energias_acumuladas_elementos_nuevo[elemento]=energias_acumuladas_elementos[elemento] + f_acumulada_elementos[elemento]'*inc_u_elemento + inc_energia_elemento
        end
        # Obtener el incremento diferencial de energia total de la zona
        inc_energia=1/2*inc_u'*inc_f    


        ### OBTENER LAS DEFORMACIONES EN LOS PUNTOS DE INTEGRACION DE CADA ELEMENTO 
        # son necesarias las deformaciones y no los incrementos para entrar a la curva tension deformacion del material
        #eps_x_ptos_int_nuevo, eps_xy_ptos_int_nuevo, eps_xz_ptos_int_nuevo = deformaciones_ptos_int(inc_u_dict,eps_x_ptos_int_acumulado, eps_xy_ptos_int_acumulado, eps_xz_ptos_int_acumulado,coords_zona_metam,conectividad_zona,radio)
        

        ### OBTENER LOS INVARIANTES Y LA DEFORMACION EQUIVALENTE DE CADA ELEMENTO
        #invariante_def=invariantes_deformaciones(tipo_invariante, eps_x_ptos_int_nuevo, eps_xy_ptos_int_nuevo, eps_xz_ptos_int_nuevo,conectividad_zona)       
        # for elemento in conectividad_zona[:,3]
        #     println(deform_equiv[elemento])
        # end
        
        ### ACTUALIZAR LOS MODULOS DE YOUNG Y OBTENER LAS ENERGIAS POR UNIDAD DE VOLUMEN DE LOS ELEMENTOS
        Youngs_tangentes_barras_actual=curva_tension_deformacion_tabla_Osgood(invariante_def_nuevo, tipo_invariante ,tensiones_pandeo)
        #energias_unitarias_elementos

        ### OBTENER LA ENERGIA TOTAL DE LA ZONA
        #energias_elementos=energia_interna_total_curva_tension_deformacion(energias_unitarias_elementos, conectividad_zona, coords_zona_metam, radio)
        #energia=sum(values(energias_elementos))


        
        # Representar las diferencias entre las energias anteriores
        println("iteracion metamaterial: ",iter," deform. equiv. max del metamaterial= ",maximum(abs.(values(invariante_def_nuevo)))," inc. energ. metamaterial= ", inc_energia) 
        
        #$$$$$
        # for elemento in conectividad_zona[:,3]
        #     println("invariante_1= ",inv_1[elemento]," energ. grafica= ",energias_elementos[elemento]," energ. lineal= ",energias_elementos_lineal[elemento])
        #     display(plots_elementos[elemento])
        #     readline()
        # end     
        
        #println();println(inv_1);println(Youngs_tangentes_barras)
    


        ### CONDICION DE FIN DE ITERACIONES 
        if abs(inc_energia-inc_energia_antiguo)/abs(inc_energia) < tolerancia1; 
            println("salida por condicion 1");  break
        else        
            if (abs(inc_energia-inc_energia_antiguo_2)/abs(inc_energia) < tolerancia1) && (abs(inc_energia-inc_energia_antiguo)/abs(inc_energia) < tolerancia2)
                println("salida por condicion 2"); break
            else
                inc_energia_antiguo_2=inc_energia_antiguo
                inc_energia_antiguo=inc_energia
            end
        end

        
    end
    #display(plot(energia_plot)) 
    if n==iter_max; println("se ha alcanzado el maximo de iteraciones en el metamaterial"); else; end

    for i in 1:length(conectividad_zona[:,3])
        elemento=conectividad_zona[i,3]
        nodo_1=conectividad_zona[i,1]; indice_1=findall(in(nodo_1),coords_zona_metam[:,4]); coords_nodo_1=coords_zona_metam[indice_1,1:3]
        nodo_2=conectividad_zona[i,2]; indice_2=findall(in(nodo_2),coords_zona_metam[:,4]); coords_nodo_2=coords_zona_metam[indice_2,1:3]
        inc_u_elemento=vcat(inc_u_dict[nodo_1],inc_u_dict[nodo_2])
        K_elemento = matriz_K_elemento(coords_nodo_1, coords_nodo_2, Youngs_tangentes_barras[elemento], radio)

        inc_f_elemento=K_elemento*inc_u_elemento
        inc_energia_elemento=1/2*inc_u_elemento'*inc_f_elemento
        inc_energias_elementos[elemento]=inc_energia_elemento

        L=norm(coords_nodo_1-coords_nodo_2)
        #invariante_def_nuevo[elemento]=invariante_def[elemento]+sqrt(inc_energia_elemento/(1/2*Youngs_tangentes_barras[elemento]*pi*radio^2*L))
        #println(elemento," ",invariante_def[elemento]," ",invariante_def_nuevo[elemento], " ", inc_energias_elementos[elemento])
        f_acumulada_elementos_nuevo[elemento]=f_acumulada_elementos[elemento] + inc_f_elemento
        energias_acumuladas_elementos_nuevo[elemento]=energias_acumuladas_elementos[elemento] + f_acumulada_elementos[elemento]'*inc_u_elemento + inc_energia_elemento
    end
    
    return inc_u, inc_f, f_acumulada_elementos_nuevo, inc_energia, inc_energias_elementos, energias_acumuladas_elementos_nuevo, Youngs_tangentes_barras, eps_x_ptos_int_nuevo, eps_xy_ptos_int_nuevo, eps_xz_ptos_int_nuevo, invariante_def_nuevo

end #function resolucion_iterativa_metamaterial



function traductores(conectividad_zona)

    N_elementos=length(conectividad_zona[:,1])
    #obtenermos los nodos a partir de la matriz de conectividad_zona    
    Nodos=vec(conectividad_zona[:,1:2]);  Nodos=sort(unique(Nodos)) #numeracion de los nodos a estudiar en el metamaterial (sin repetir y ordenados)
    N_gdl_zona=6*length(Nodos)
    gdl_traductor_zona=zeros(Int64,N_gdl_zona,2); elementos_traductor_zona=zeros(Int64, N_elementos, 2)

    i=1
    for nodo in Nodos 
        gdl_traductor_zona[6*i-5:6*i,1]=collect(6*nodo-5:6*nodo) # gdl a estudiar en la region (6 por nodo)
        gdl_traductor_zona[6*i-5:6*i,2]=collect(6*i-5:6*i) # esto sirve para asociar los gdl a estudiar con una sucesion de numeros ordenada entre 1 y el numero maximo de gdl
                                            # ...actúa como de traductor
        i=i+1
    end


    for elemento in 1:N_elementos
        elementos_traductor_zona[elemento,1]=conectividad_zona[elemento,3]
        elementos_traductor_zona[elemento,2]=elemento
    end

    return gdl_traductor_zona, elementos_traductor_zona
    
end # function traductores


function generar_diccionario_u(inc_u,coords_zona_metam) # permite obtener el nodo del total y su correspondiente vector de desplazamientos en el nodo
    N_nodos=length(coords_zona_metam[:,1])
    inc_u_dict=Dict()

    for i in 1:N_nodos
        # obtener el valor del nodo global
        Nodo_metam_global=coords_zona_metam[i,4]

        inc_u_dict[Nodo_metam_global]=inc_u[6*i-5:6*i]
    end

    return inc_u_dict

end # function generar_diccionario_u


function generar_diccionarios(inc_u,inc_f, coords_zona_metam) # permite obtener el nodo del total y su correspondiente vector de desplazamientos en el nodo
    N_nodos=length(coords_zona_metam[:,1])
    inc_desplaz_dict=Dict(); inc_giros_dict=Dict(); inc_fuerzas_dict=Dict(); inc_momentos_dict=Dict()

    for i in 1:N_nodos
        # obtener el valor del nodo global
        Nodo_metam_global=coords_zona_metam[i,4]

        inc_desplaz_dict[Nodo_metam_global]=inc_u[6*i-5:6*i-3]
        inc_giros_dict[Nodo_metam_global]=inc_u[6*i-2:6*i]
        inc_fuerzas_dict[Nodo_metam_global]=inc_f[6*i-5:6*i-3]
        inc_momentos_dict[Nodo_metam_global]=inc_f[6*i-2:6*i]

    end

    return inc_desplaz_dict, inc_giros_dict, inc_fuerzas_dict, inc_momentos_dict

end # function generar_diccionario_u


function tensiones_criticas_pandeo(coords_zona_metam,conectividad_zona_metam,young_material, radio)

    alpha=0.5 # parametro para determinar como modelamos el pandeo --> entre 0 y 1
              # (alpha) multiplica al valor de pandeo de una viga biempotrada
              # (1-alpha) multiplica al valor de pandeo de una viga simplemente apoyada

    sigmas_pandeo=Dict()
    N_elementos=length(conectividad_zona_metam[:,3])

    for i in 1:N_elementos

        elemento=conectividad_zona_metam[i,3]

        indice_1=findall(in(conectividad_zona_metam[i,1]), coords_zona_metam[:,4])[1]
        coords_nodo_1= coords_zona_metam[indice_1,1:3]
        indice_2=findall(in(conectividad_zona_metam[i,2]), coords_zona_metam[:,4])[1]
        coords_nodo_2= coords_zona_metam[indice_2,1:3]

        longitud=norm(coords_nodo_2-coords_nodo_1)
        Iy= 1/4*pi*(radio^2)^2  #momento de inercia respecto al eje Y. Al ser circular da igual el Y que el Z
        Area=pi*radio^2

        sigma_pandeo_superior=pi^2*young_material*Iy/(0.5*longitud^2*Area) # Modelamos el elemento como una viga biempotrada solo trabajando a compresion (es la fuerza principal)
        sigma_pandeo_inferior=pi^2*young_material*Iy/(longitud^2*Area)  # modelamos como una viga simplemente apoyada
        
        # la tension de pandeo estará entre las dos anteriores
        sigmas_pandeo[elemento]=alpha*sigma_pandeo_superior+(1-alpha)*sigma_pandeo_inferior

    end
    #$$$
    #println("tensiones_pandeo");println(sigmas_pandeo)

    return sigmas_pandeo

end


function error_fuerzas(inc_f, gdl_zona_restricciones_no_fijas, gdl_zona_metam)
    
    N_restricciones_no_fijas=length(gdl_zona_restricciones_no_fijas)
    f_no_fijo=zeros(N_restricciones_no_fijas)
    

    for i in 1:N_restricciones_no_fijas
        gdl_restring_no_fijo= gdl_zona_restricciones_no_fijas[i]

        # obtener la numeracion local de este gdl
        indice_local=findall(in(gdl_restring_no_fijo),gdl_zona_metam)[1]

        # obtener la fuerza asociada a ese nodo
        f_no_fijo[i]=inc_f[indice_local]
    end

    return f_no_fijo

end


function derivadas_nodos_no_fijos(f_no_fijo_global, coords_metam_hexaedros_exacto, conectiv_metam_hexaedros_exacto, gdl_restricciones_no_fijas_global, gdl_zona_restricciones_no_fijas_hexaedros, gdl_zona_metam_hexaedros, gdl_zona_libres_hexaedros, Youngs_tangentes_memoria_metam, radio, longitud_caracteristica)
    
    N_restricciones_no_fijas=length(gdl_restricciones_no_fijas_global)
    df_du_no_fijo_global=zeros(N_restricciones_no_fijas)
    df_du_no_fijo_zona_matriz=Dict()
    

    for i in 1:N_restricciones_no_fijas
        gdl_no_fijo_u=gdl_restricciones_no_fijas_global[i]
        indice_u_global=i

        for elemento_solido in 1:length(coords_metam_hexaedros_exacto)
            gdl_zona_restricciones_no_fijas=gdl_zona_restricciones_no_fijas_hexaedros[elemento_solido]
            indice_u_local=findall(in(gdl_no_fijo_u), gdl_zona_restricciones_no_fijas)

            # Hay que ver si el gdl que estudiamos se encuentra en el elemento o no
            if indice_u_local != []
                # si ya hemos calculado la zona, no volvemos a calcular la inversa
                N_restricciones_no_fijas_zona=length(gdl_zona_restricciones_no_fijas)
                
                if haskey(df_du_no_fijo_zona_matriz, elemento_solido)==true
                else
                    K_zona_metam =matriz_K_estructura(conectiv_metam_hexaedros_exacto[elemento_solido], coords_metam_hexaedros_exacto[elemento_solido], Youngs_tangentes_memoria_metam[elemento_solido], radio)
                    gdl_zona_libres=gdl_zona_libres_hexaedros[elemento_solido]            
                    gdl_zona_metam=gdl_zona_metam_hexaedros[elemento_solido]
                    gdl_restric_no_fijas_local=indexin(gdl_zona_restricciones_no_fijas, gdl_zona_metam)
                    gdl_libres_local=indexin(gdl_zona_libres, gdl_zona_metam)

                    # obtener las derivadas de las fuerzas en estos nodos con respecto al desplazamiento de esos nodos    
                    df_du_no_fijo_zona_matriz[elemento_solido]=K_zona_metam[gdl_restric_no_fijas_local, gdl_restric_no_fijas_local]-K_zona_metam[gdl_restric_no_fijas_local,gdl_libres_local]*(inv(K_zona_metam[gdl_libres_local,gdl_libres_local])*K_zona_metam[gdl_libres_local,gdl_restric_no_fijas_local])
                    # for i in 1:N_restricciones_no_fijas_zona
                    #     println(gdl_zona_restricciones_no_fijas[i])
                    #     println(elemento_solido)
                    #     println(df_du_no_fijo_zona_matriz[elemento_solido][:,i])
                    #     println()
                    # end
                    # readline()
                end           
                indice_u_local=indice_u_local[1]

                # Factor para tener en cuenta el diferente orden de magnitud de giros y desplazamientos
                k=gdl_no_fijo_u%6
                if k!=0 && k<=3
                    factor_u=1.0
                else
                    factor_u=1/longitud_caracteristica 
                end

                for j in 1:N_restricciones_no_fijas_zona
                    gdl_no_fijo_f=gdl_zona_restricciones_no_fijas[j]
                    indice_f_global=findall(in(gdl_no_fijo_f),gdl_restricciones_no_fijas_global)[1]

                    # Factor para tener en cuenta el diferente orden de magnitud de fueras y momentos
                    k=gdl_no_fijo_f%6
                    if k!=0 && k<=3
                        factor_f=1.0
                    else
                        factor_f=1/longitud_caracteristica
                    end

                    df_du_no_fijo_global[indice_u_global] = df_du_no_fijo_global[indice_u_global] + factor_u*factor_f*2*f_no_fijo_global[indice_f_global]*df_du_no_fijo_zona_matriz[elemento_solido][j,indice_u_local]
                
                    #println(gdl_no_fijo_u," ",gdl_no_fijo_f," ", f_no_fijo_global[indice_f_global]," ",df_du_no_fijo_zona_matriz[elemento_solido][j,indice_u_local]," ",df_du_no_fijo_global[indice_u_global]);readline()
                end
            else
            end
        end
    end

    return df_du_no_fijo_global

end


# igual que la funcion anterior pero indicando el hexaedro al que le corresponde la zona estudiada 
function derivadas_nodos_no_fijos_elemento(elemento, df_du_no_fijo_zona_matriz, f_no_fijo_global, hexaedros_zona, coords_metam_hexaedros_exacto, conectiv_metam_hexaedros_exacto, gdl_restricciones_no_fijas_global, gdl_zona_restricciones_no_fijas_hexaedros, gdl_zona_metam_hexaedros, gdl_zona_libres_hexaedros, Youngs_tangentes_memoria_metam, radio, longitud_caracteristica)
    
    N_restricciones_no_fijas=length(gdl_restricciones_no_fijas_global)
    df_du_no_fijo_global=zeros(N_restricciones_no_fijas)
    #df_du_no_fijo_zona_matriz=Dict()
    

    for i in 1:N_restricciones_no_fijas
        gdl_no_fijo_u=gdl_restricciones_no_fijas_global[i]
        indice_u_global=i

        for elemento_solido in hexaedros_zona
            gdl_zona_restricciones_no_fijas=gdl_zona_restricciones_no_fijas_hexaedros[elemento, elemento_solido]
            indice_u_local=findall(in(gdl_no_fijo_u), gdl_zona_restricciones_no_fijas)

            # Hay que ver si el gdl que estudiamos se encuentra en el elemento o no
            if indice_u_local != []
                # si ya hemos calculado la zona, no volvemos a calcular la inversa
                N_restricciones_no_fijas_zona=length(gdl_zona_restricciones_no_fijas)
                
                if haskey(df_du_no_fijo_zona_matriz, elemento_solido)==true                 
                else
                    K_zona_metam =matriz_K_estructura(conectiv_metam_hexaedros_exacto[elemento_solido], coords_metam_hexaedros_exacto[elemento_solido], Youngs_tangentes_memoria_metam[elemento, elemento_solido], radio)
                    gdl_zona_libres=gdl_zona_libres_hexaedros[elemento, elemento_solido]            
                    gdl_zona_metam=gdl_zona_metam_hexaedros[elemento_solido]
                    gdl_restric_no_fijas_local=indexin(gdl_zona_restricciones_no_fijas, gdl_zona_metam)
                    gdl_libres_local=indexin(gdl_zona_libres, gdl_zona_metam)

                    # obtener las derivadas de las fuerzas en estos nodos con respecto al desplazamiento de esos nodos    
                    df_du_no_fijo_zona_matriz[elemento_solido]=K_zona_metam[gdl_restric_no_fijas_local, gdl_restric_no_fijas_local]-K_zona_metam[gdl_restric_no_fijas_local,gdl_libres_local]*(inv(K_zona_metam[gdl_libres_local,gdl_libres_local])*K_zona_metam[gdl_libres_local,gdl_restric_no_fijas_local])
                    # for i in 1:N_restricciones_no_fijas_zona
                    #     println(gdl_zona_restricciones_no_fijas[i])
                    #     println(elemento_solido)
                    #     println(df_du_no_fijo_zona_matriz[elemento_solido][:,i])
                    #     println()
                    # end
                    # readline()
                end           
                indice_u_local=indice_u_local[1]

                # Factor para tener en cuenta el diferente orden de magnitud de giros y desplazamientos
                k=gdl_no_fijo_u%6
                if k!=0 && k<=3
                    factor_u=1.0
                else
                    factor_u=1/longitud_caracteristica
                end

                for j in 1:N_restricciones_no_fijas_zona
                    gdl_no_fijo_f=gdl_zona_restricciones_no_fijas[j]
                    indice_f_global=findall(in(gdl_no_fijo_f),gdl_restricciones_no_fijas_global)[1]

                    # Factor para tener en cuenta el diferente orden de magnitud de fueras y momentos
                    k=gdl_no_fijo_f%6
                    if k!=0 && k<=3
                        factor_f=1.0
                    else
                        factor_f=1/longitud_caracteristica
                    end

                    df_du_no_fijo_global[indice_u_global] = df_du_no_fijo_global[indice_u_global] + factor_u*factor_f*2*f_no_fijo_global[indice_f_global]*df_du_no_fijo_zona_matriz[elemento_solido][j,indice_u_local]
                
                    #println(gdl_no_fijo_u," ",gdl_no_fijo_f," ", f_no_fijo_global[indice_f_global]," ",df_du_no_fijo_zona_matriz[elemento_solido][j,indice_u_local]," ",df_du_no_fijo_global[indice_u_global]);readline()
                end
            else
            end
        end
    end

    return df_du_no_fijo_global, df_du_no_fijo_zona_matriz

end


function derivadas_nodos_no_fijos_K(coords_metam_hexaedros_exacto, conectiv_metam_hexaedros_exacto, gdl_restricciones_no_fijas_global, gdl_no_fijos_zona_hexaedros, gdl_restring_no_fijos_hexaedro_local, Youngs_tangentes_memoria_metam, radio, K_step, gradiente_fuerzas_nodos, inc_valor_restricciones_no_fijas_total)
    
    N_restricciones_no_fijas=length(gdl_restricciones_no_fijas_global)
    derivadas_u=Dict()


    for elemento_solido in 1:length(coords_metam_hexaedros_exacto)
        println(elemento_solido)
        K_zona_metam =matriz_K_estructura(conectiv_metam_hexaedros_exacto[elemento_solido], coords_metam_hexaedros_exacto[elemento_solido], Youngs_tangentes_memoria_metam[elemento_solido], radio)
        indices=indexin(gdl_no_fijos_zona_hexaedros[elemento_solido], gdl_restricciones_no_fijas_global)

        for j in 1:length(gdl_no_fijos_zona_hexaedros[elemento_solido])
            gdl_no_fijo_global=gdl_no_fijos_zona_hexaedros[elemento_solido][j]
            gdl_no_fijo_local=gdl_restring_no_fijos_hexaedro_local[elemento_solido][j]

            if haskey(derivadas_u, gdl_no_fijo_global)==false
                derivadas_u[gdl_no_fijo_global]=zeros(N_restricciones_no_fijas)
                derivadas_u[gdl_no_fijo_global][indices]=K_zona_metam[gdl_restring_no_fijos_hexaedro_local[elemento_solido], gdl_no_fijo_local]
            else
                derivadas_u[gdl_no_fijo_global][indices]=derivadas_u[gdl_no_fijo_global][indices] + K_zona_metam[gdl_restring_no_fijos_hexaedro_local[elemento_solido], gdl_no_fijo_local]
            end
        end
    end

    
    for j in 1:N_restricciones_no_fijas
        gdl_no_fijo=gdl_restricciones_no_fijas_global[j]
        inc_valor_restricciones_no_fijas_total[j]=inc_valor_restricciones_no_fijas_total[j]-K_step[j]*derivadas_u[gdl_no_fijo]'*gradiente_fuerzas_nodos
    end

    return inc_valor_restricciones_no_fijas_total

end


function energia_barras_aisladas(coords_nodo_1, coords_nodo_2, u_elemento, tension_pandeo, iter_desplaz, radio)

    N_incrementos_u=iter_desplaz   # Para poder estudiar bien la parte no lineal 
    inc_u_elemento=u_elemento/N_incrementos_u
    poisson=Poisson(); Young_elemento=Young_material()
    energia_barra_acumulada=0; f_barra_acumulada=zeros(12)
    eps_x_ptos_int_acumulado=[0,0]; eps_xy_ptos_int_acumulado=[0,0]; eps_xz_ptos_int_acumulado=[0,0];
    inc_energia_barra=0

    K_elemento_por_unidad_de_young=matriz_K_elemento(coords_nodo_1, coords_nodo_2, Young_elemento, radio)/Young_elemento

    ######################## INCREMENTO DE ENERGIA DEL ELEMENTO POR UNIDAD DE YOUNG
    inc_energia_barra_por_unidad_de_young=1/2*inc_u_elemento'*K_elemento_por_unidad_de_young*inc_u_elemento


    for i in 1:N_incrementos_u
        
        ################################### DEFORMACIONES DE LA BARRA

        chi_ptos_int=[-1/sqrt(3),1/sqrt(3)]   #coordenada naturales de los puntos de integración

        R_12_12=ejes_locales_a_globales(coords_nodo_1,coords_nodo_2) #matriz para pasar de coordenadas locales a globales

        u_elemento_loc=R_12_12'*u_elemento
        u_x1=u_elemento_loc[1];             u_x2=u_elemento_loc[7];
        u_y1=u_elemento_loc[2];             u_y2=u_elemento_loc[8];
        u_z1=u_elemento_loc[3];             u_z2=u_elemento_loc[9];
        theta_x1=u_elemento_loc[4];         theta_x2=u_elemento_loc[10];
        theta_y1=u_elemento_loc[5];         theta_y2=u_elemento_loc[11];
        theta_z1=u_elemento_loc[6];         theta_z2=u_elemento_loc[12];


        # Jacobiano del elemento
        L = norm(coords_nodo_2-coords_nodo_1)  # longitud del elemento
        Jac=L/2; inv_Jac=1/Jac

        y=radio; z=radio; r=radio             #para obtener el valor máximo de deformacion en los puntos de integracion

        #$$$$
        inc_eps_x_axil=0;

        eps_x_ptos_int_acumulad=eps_x_ptos_int_acumulado; eps_xy_ptos_int_acumulad=eps_xy_ptos_int_acumulado; eps_xz_ptos_int_acumulad=eps_xz_ptos_int_acumulado
        eps_x_ptos_int=[]; eps_xy_ptos_int=[]; eps_xz_ptos_int=[]

        for i in 1:length(chi_ptos_int)

            ## INCREMENTO DEFORMACION epsilon_x POR PARTE AXIL (2 puntos de integracion)--> salvo que sea nulo, el axil suele ser el dominante
            N=1/2*[1-chi_ptos_int[i],1+chi_ptos_int[i]]
            B=inv_Jac*1/2*[-1,1]
            inc_eps_x_axil=dot(B,[u_x1,u_x2])

            ### INCREMENTO DEFORMACION epsilon_x POR PARTE DE FLEXION EN Y y en Z LOCAL a una distancia del eje (radio) (viga hermitica)
            N_herm_1=-0.75*chi_ptos_int[i]+0.25*chi_ptos_int[i]^3+0.5
            N_herm_2=-0.25*chi_ptos_int[i]-0.25*chi_ptos_int[i]^2+0.25*chi_ptos_int[i]^3+0.25
            N_herm_3=0.75*chi_ptos_int[i]-0.25*chi_ptos_int[i]^3+0.5
            N_herm_4=-0.25*chi_ptos_int[i]+0.25*chi_ptos_int[i]^2+0.25*chi_ptos_int[i]^3-0.25
            N_herm=[N_herm_1,Jac*N_herm_2,N_herm_3,Jac*N_herm_4]
            B_herm=inv_Jac^2*[6/4*chi_ptos_int[i], Jac*(-2+6*chi_ptos_int[i])/4, -6/4*chi_ptos_int[i], Jac*(2+6*chi_ptos_int[i])/4]

            # Hay que determinar el punto de la seccion de los ptos de integracion (coordenadas Y y Z) donde la deformacion es maxima (será el mas critico)
            inc_eps_x_flex_y_aux=dot(B_herm,[u_y1,theta_z1,u_y2,theta_z2])
            inc_eps_x_flex_z_aux=dot(B_herm,[u_z1,theta_y1,u_z2,theta_y2])

            
            alpha=atan(-inc_eps_x_flex_z_aux,inc_eps_x_flex_y_aux) #resolviendo: maximo(-radio*cos(alpha)*inc_eps_x_flex_y_aux+radio*sin(alpha)*inc_eps_x_flex_z_aux)

            inc_eps_x_flex_y=radio*cos(alpha)*inc_eps_x_flex_y_aux
            inc_eps_x_flex_z=radio*sin(alpha)*inc_eps_x_flex_z_aux
            inc_eps_x_flex=-inc_eps_x_flex_y + inc_eps_x_flex_z


            ### Se debe tener en cuenta si eps_axil acumulado en el elemento es positivo o negativo a la hora de seleccionar el inc_epsilon de flexion que hace mas critica a la viga
            if inc_eps_x_axil > 0
                push!(eps_x_ptos_int, eps_x_ptos_int_acumulad[i]+inc_eps_x_axil+abs(inc_eps_x_flex))
            elseif inc_eps_x_axil < 0
                push!(eps_x_ptos_int, eps_x_ptos_int_acumulad[i]+inc_eps_x_axil-abs(inc_eps_x_flex))
            else #si el axil es cero se coge el positivo (de traccion) (criterio arbitrario)
                push!(eps_x_ptos_int, eps_x_ptos_int_acumulad[i]+abs(inc_eps_x_flex))
            end


            ### INCREMENTOS DE DEFORMACION gamma_x_r a una distancia del eje (r) POR EL GIRO ENTORNO AL EJE X
            inc_gamma_x_r=r*dot(B,[theta_x1,theta_x2])

            ### a partir del gamma_x_r, obtenemos gamma_xy y gamma_xz en el punto donde es máximo el inc_epsilon_x --> eso se determina por el alpha
            inc_gamma_xy=-inc_gamma_x_r*sin(alpha) #el signo no importa porque se va a usar para los invariantes y alli se usa al cuadrado
            inc_gamma_xz=inc_gamma_x_r*cos(alpha)
            push!(eps_xy_ptos_int, eps_xy_ptos_int_acumulad[i]+1/2*inc_gamma_xy)
            push!(eps_xz_ptos_int, eps_xz_ptos_int_acumulad[i]+1/2*inc_gamma_xz)
        
        end


        ######################### CALCULO DE LA DEFORMACION EQUIVALENTE

        invariantes_barra=[]
        chi_min=-sqrt(3); chi_max=sqrt(3);
        for chi in range(chi_min,chi_max,length=5)  # se calculan los invariantes en varios puntos para luego coger el que dé mayores valores; mas puntos = mas precision para coger el maximo 
            N=1/2*[1-chi,1+chi]
            eps_x=dot(N,eps_x_ptos_int)
            eps_xy=dot(N,eps_xy_ptos_int)
            eps_xz=dot(N,eps_xz_ptos_int)
            eps_y=-poisson*eps_x; eps_z=-poisson*eps_x #estimamos el eps_y y el eps_x

            signo=sign(eps_x)

            deformacion_equivalente=sqrt(1/((1+poisson)*(1-2*poisson))*((1-poisson)*(eps_x^2+eps_y^2+eps_z^2)+2*poisson*(eps_x*eps_y + eps_x*eps_z + eps_z*eps_y))+4/(2*(1+poisson))*(eps_xy^2+eps_xz^2))
            append!(invariantes_barra, signo*deformacion_equivalente)
        end

        invariante_max=maximum(invariantes_barra);  invariante_min=minimum(invariantes_barra)

        if abs(invariante_max) >= abs(invariante_min);   deform_equiv=invariante_max
        else;                                            deform_equiv=invariante_min
        end


        ######################## ACTUALIZAR EL YOUNG TANGENTE
        Youngs_tangentes_barras, energias_unitarias_elementos= curva_tension_deformacion_tabla_Osgood_1_elemento(deform_equiv, "deform_equiv", tension_pandeo)


        ######################## INCREMENTO DE ENERGIA Y FUERZA DEL ELEMENTO
        inc_energia_barra= inc_energia_barra_por_unidad_de_young * Young_elemento
        energia_barra_acumulada= energia_barra_acumulada + inc_energia_barra + f_barra_acumulada'*inc_u_elemento
        f_barra_acumulada= f_barra_acumulada + K_elemento_por_unidad_de_young*Young_elemento*inc_u_elemento

    end

    return energia_barra_acumulada, inc_energia_barra


end


end # module