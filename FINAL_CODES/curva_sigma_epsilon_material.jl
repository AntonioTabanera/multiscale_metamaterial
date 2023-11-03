# MODULO PARA IMPLEMENTAR LA CURVA DE TENSION-DEFORMACION DEL MATERIAL BASE EN CADA INSTANTE
module curva_sigma_epsilon_material

using DelimitedFiles

export curva_tension_deformacion_Osgood_analitica, generar_tabla_tension_deformacion_Osgood, curva_tension_deformacion_tabla_Osgood, curva_tension_deformacion_tabla_Osgood_1_elemento, Poisson, Young_material


function curva_tension_deformacion_Osgood_analitica(epsilons_ingenieril, nombre_invariante, tensiones_pandeo) # MODULO SECANTE --> obtener el modulo elastico equivalente (dañado) 

    energias_especificas=Dict(); Youngs_tangentes_ing=Dict() # Modulo de Young tangente de la curva ingenieril
    poisson=Poisson()

    # CON ESTE MODELO SE USAN DEFORMACIONES REALES (EN TODOS LOS CALCULOS HEMOS USADO LAS INGENIERILES) 
    # las curvas tension real- deformacion real son practicamente iguales en traccion que en compresion

    ### PARAMETROS DEL MATERIAL (acero uns30403 rhs)
    Young=194000  # MPa
    sigma_02_real=445 # MPa
    sigma_ultimo_real=730 # MPa
    epsilon_ultimo_real=0.51
    n=4.85
    epsilon_fractura=0.53  # valor inventado pero tipico en aceros frente al epsilon_ultimo_real
    sigma_real_limite_compresion=1000
 
    ### PARAMETROS DERIVADOS
    epsilon_02_real=sigma_02_real/Young+0.002
    e=sigma_02_real/Young
    m=1+3.5*sigma_02_real/sigma_ultimo_real
    Young_02=Young/(1+0.002*n/e) # modulo de young tangente (la derivada) en 0.2      

    r=Young*epsilon_02_real/sigma_02_real
    r_2=Young_02*epsilon_02_real/sigma_02_real
    r_asterisco=Young_02*(epsilon_ultimo_real-epsilon_02_real)/(sigma_ultimo_real-sigma_02_real)
    Young_u=Young_02/(1+(r_asterisco-1)*m)
    r_u=Young_u*(epsilon_ultimo_real-epsilon_02_real)/(sigma_ultimo_real-sigma_02_real)    
    p=r*(1-r_2)/(r-1)
    p_asterisco=r_asterisco*(1-r_u)/(r_asterisco-1)
    epsilon_ult_p_asterisco=epsilon_ultimo_real-epsilon_02_real-(sigma_ultimo_real-sigma_02_real)/Young_02
    epsilon_ult_normalizado=epsilon_ultimo_real/epsilon_02_real
    epsilon_ult_ingenieril=exp(epsilon_ultimo_real)-1
    epsilon_fractura_ingenieril=exp(epsilon_fractura)-1


    ### MODELO APROXIMADO A PARTIR DE LA INVERSION DEL MODELO DE RAMBERG OSGOOD
    N_elementos=length(epsilons_ingenieril)
    elementos=collect(keys(epsilons_ingenieril))
    

    for i in 1:N_elementos

        # El dato es uno de los invariantes --> obtendremos el epsilon_x del ensayo de traccion para el cual se obtendría el mismo...
        #... valor de invariante --> con ese epsilon_x ingenieril pasaremos al epsilon_x real y de ahí entraremos en el modelo analitico de Ramberg-Osgood

        if nombre_invariante=="inv_1"; epsilon_ingenieril=epsilons_ingenieril[elementos[i]]/(1-2*poisson)
        elseif nombre_invariante=="inv_2"; epsilon_ingenieril=epsilons_ingenieril[elementos[i]]/sqrt(2*poisson-poisson^2)
        elseif nombre_invariante=="inv_3"; epsilon_ingenieril=epsilons_ingenieril[elementos[i]]/cbrt(poisson^2)
        elseif nombre_invariante=="deform_equiv"; epsilon_ingenieril=epsilons_ingenieril[elementos[i]]
        elseif nombre_invariante=="deform_equiv_2"; epsilon_ingenieril=epsilons_ingenieril[elementos[i]]
        else println("no se ha introducido un criterio para entrar en la curva tension-deformacion")
        end

        
        # DEFORMACION_NULA
        if epsilon_ingenieril==0; Youngs_tangentes_ing[elementos[i]]=Young; energias_especificas[elementos[i]]=0.0

        # TRACCION:   
        elseif epsilon_ingenieril > 0

            # PARTE ELASTICA + PLASTICA  
            if epsilon_ingenieril <= epsilon_ult_ingenieril  
                epsilon_real=log(1+epsilon_ingenieril)
                epsilon_real_normalizado=epsilon_real/epsilon_02_real
                
                if epsilon_real_normalizado <= 1
                    sigma_real=(r*epsilon_real_normalizado/(1+(r-1)*epsilon_real_normalizado^p))*sigma_02_real
                    energias_especificas[elementos[i]]=sigma_real*epsilon_real-(sigma_real^2/(2*Young)+0.002/(n+1)*sigma_real^(n+1)/sigma_02_real^n)
                    Youngs_tangentes_ing[elementos[i]]=(1/(1/Young+0.002*n*(sigma_real^(n-1))/(sigma_02_real^n))-sigma_real)*1/((1+epsilon_ingenieril)^2)
                else
                    sigma_real=(1+r_2*(epsilon_real_normalizado-1)/(1+(r_asterisco-1)*((epsilon_real_normalizado-1)/(epsilon_ult_normalizado-1))^p_asterisco))*sigma_02_real
                    energias_especificas[elementos[i]]=sigma_real*epsilon_real-(epsilon_02_real*sigma_real+sigma_real^2/(2*Young_02)-sigma_02_real*sigma_real/Young_02+epsilon_ult_p_asterisco/(m+1)*(sigma_real-sigma_02_real)^(m+1)/(sigma_ultimo_real-sigma_02_real)^m-epsilon_02_real*sigma_02_real+0.002/(n+1)*sigma_02_real+sigma_02_real^2/2*(1/Young+1/Young_02))
                    Youngs_tangentes_ing[elementos[i]]=(1/(1/Young_02+epsilon_ult_p_asterisco*m*((sigma_real-sigma_02_real)^(m-1))/((sigma_ultimo_real-sigma_02_real)^m))-sigma_real)/((1+epsilon_ingenieril)^2)
                end       
                sigma_ingenieril=sigma_real/(1+epsilon_ingenieril)


            # PARTE DE FRACTURA (ABLANDAMIENTO) -- modelado como una linea recta
            elseif (epsilon_ingenieril > epsilon_ult_ingenieril) && (epsilon_ingenieril <= epsilon_fractura_ingenieril)
                energia_acumulada=sigma_ultimo_real*epsilon_ultimo_real-(epsilon_02_real*sigma_ultimo_real+sigma_ultimo_real^2/(2*Young_02)-sigma_02_real*sigma_ultimo_real/Young_02+epsilon_ult_p_asterisco/(m+1)*(sigma_ultimo_real-sigma_02_real)-0.002/(n+1)*sigma_02_real)
                sigma_ingenieril_ult=sigma_ultimo_real/(1+epsilon_ult_ingenieril)
                sigma_ingenieril=sigma_ingenieril_ult-(epsilon_ingenieril-epsilon_ult_ingenieril)/(epsilon_fractura_ingenieril-epsilon_ult_ingenieril)*sigma_ingenieril_ult
                Youngs_tangentes_ing[elementos[i]]=-(epsilon_ingenieril-epsilon_ult_ingenieril)/(epsilon_fractura_ingenieril-epsilon_ult_ingenieril)
                energias_especificas[elementos[i]]=energia_acumulada + sigma_ingenieril_ult*epsilon_ingenieril-(1/2*epsilon_ingenieril^2-epsilon_ult_ingenieril*epsilon_ingenieril)/(epsilon_fractura_ingenieril-epsilon_ult_ingenieril)*sigma_ingenieril_ult

            # ELEMENTO ROTO
            elseif (epsilon_ingenieril > epsilon_fractura)
                Youngs_tangentes_ing[elementos[i]]=0.001 #valor muy pequeño de modulo de Young pero que no sea cero para evitar fallos al invertir las matrices
                energias_especificas[elementos[i]]=0

            else
            end

        # COMPRESION: (a compresion suponemos que no rompen. O pandea o aguanta)
        else
            epsilon_ingenieril=abs(epsilon_ingenieril)
            sigma_ing_pandeo=tensiones_pandeo[elementos[i]]
        

            # PARTE ELASTICA + PLASTICA + PANDEO
            epsilon_real=log(1+epsilon_ingenieril)
            epsilon_real_normalizado=epsilon_real/epsilon_02_real
            
            # Primer tramo de la curva
            if epsilon_real_normalizado <= 1
                sigma_real=(r*epsilon_real_normalizado/(1+(r-1)*epsilon_real_normalizado^p))*sigma_02_real
                sigma_ingenieril=sigma_real/(1+epsilon_ingenieril)

                # COMPROBAR PANDEO
                if sigma_ingenieril <= sigma_ing_pandeo # No hay pandeo
                    energias_especificas[elementos[i]]=sigma_real*epsilon_real-(sigma_real^2/(2*Young)+0.002/(n+1)*sigma_real^(n+1)/sigma_02_real^n)
                    Youngs_tangentes_ing[elementos[i]]=(1/(1/Young+0.002*n*(sigma_real^(n-1))/(sigma_02_real^n))-sigma_real)*1/((1+epsilon_ingenieril)^2)

                else # Hay pandeo

                    ### Determinar la deformacion critica a partir de la cual hay pandeo --> proceso iterativo
                    iter_max=30; tolerancia=1e-15; sigma_pandeo_real_antiguo=0; sigma_pandeo_real=sigma_ing_pandeo # inicializamos para el proceso iterativo
                    for k in 1:iter_max
                        sigma_pandeo_real=sigma_ing_pandeo*exp(sigma_pandeo_real/Young+0.002*(sigma_pandeo_real/sigma_02_real)^n)
                        if abs(sigma_pandeo_real-sigma_pandeo_real_antiguo)/sigma_pandeo_real < tolerancia
                            break
                        else
                            sigma_pandeo_real_antiguo=sigma_pandeo_real
                        end
                    end

                    epsilon_critico_pandeo_real=sigma_pandeo_real/Young+0.002*(sigma_pandeo_real/sigma_02_real)^n
                    energia_critica_pandeo=sigma_pandeo_real*epsilon_critico_pandeo_real-(sigma_pandeo_real^2/(2*Young)+0.002/(n+1)*sigma_pandeo_real^(n+1)/sigma_02_real^n)
       
                    sigma_ingenieril=sigma_ing_pandeo
                    energias_especificas[elementos[i]]=sigma_pandeo_real*(epsilon_real-epsilon_critico_pandeo_real)+ energia_critica_pandeo
                    Youngs_tangentes_ing[elementos[i]]=0.0001
                end

            
            else    # Si el valor de la deformacion de entrada esta en el segundo tramo
                sigma_real=(1+r_2*(epsilon_real_normalizado-1)/(1+(r_asterisco-1)*((epsilon_real_normalizado-1)/(epsilon_ult_normalizado-1))^p_asterisco))*sigma_02_real
                sigma_ingenieril=sigma_real/(1+epsilon_ingenieril)

                # COMPROBAR PANDEO
                if sigma_ingenieril <= sigma_ing_pandeo # No hay pandeo
                    energias_especificas[elementos[i]]=sigma_real*epsilon_real-(epsilon_02_real*sigma_real+sigma_real^2/(2*Young_02)-sigma_02_real*sigma_real/Young_02+epsilon_ult_p_asterisco/(m+1)*(sigma_real-sigma_02_real)^(m+1)/(sigma_ultimo_real-sigma_02_real)^m-epsilon_02_real*sigma_02_real+0.002/(n+1)*sigma_02_real+sigma_02_real^2/2*(1/Young+1/Young_02))
                    Youngs_tangentes_ing[elementos[i]]=(1/(1/Young_02+epsilon_ult_p_asterisco*m*((sigma_real-sigma_02_real)^(m-1))/((sigma_ultimo_real-sigma_02_real)^n))-sigma_real)/((1+epsilon_ingenieril)^2)

                else # Hay pandeo

                    ### Determinar la deformacion critica a partir de la cual hay pandeo --> proceso iterativo        
                    iter_max=30; tolerancia=1e-15; sigma_pandeo_real_antiguo=0; sigma_pandeo_real=sigma_ing_pandeo # inicializamos para el proceso iterativo
                    # comenzamos las iteraciones suponiendo que pandea antes de sigma_02_real
                    for k in 1:iter_max
                        sigma_pandeo_real=sigma_ing_pandeo*exp(sigma_pandeo_real/Young+0.002*(sigma_pandeo_real/sigma_02_real)^n)
                        if abs(sigma_pandeo_real-sigma_pandeo_real_antiguo)/sigma_pandeo_real < tolerancia
                            break
                        else
                            sigma_pandeo_real_antiguo=sigma_pandeo_real
                        end
                    end

                    # si el valor obtenido es menor que sigma_02_real, habremos acertado con la hipotesis de partida
                    if sigma_pandeo_real<sigma_02_real 
                        epsilon_critico_pandeo_real=sigma_pandeo_real/Young+0.002*(sigma_pandeo_real/sigma_02_real)^n
                        energia_critica_pandeo=sigma_pandeo_real*epsilon_critico_pandeo_real-(sigma_pandeo_real^2/(2*Young)+0.002/(n+1)*sigma_pandeo_real^(n+1)/sigma_02_real^n)
                    
                    else #si el valor obtenido es menor que sigma_02_real, debemos estudiar en la otra region del grafico
                        iter_max=30; tolerancia=1e-15; sigma_pandeo_real_antiguo=0; sigma_pandeo_real=sigma_ing_pandeo # inicializamos para el proceso iterativo
                        for k in 1:iter_max
                            sigma_pandeo_real=sigma_ing_pandeo*exp(epsilon_02_real+(sigma_pandeo_real-sigma_02_real)/Young_02+epsilon_ult_p_asterisco*(sigma_pandeo_real-sigma_02_real)^m/(sigma_ultimo_real-sigma_02_real)^m)
                            if abs(sigma_pandeo_real-sigma_pandeo_real_antiguo)/sigma_pandeo_real < tolerancia
                                break
                            else
                                sigma_pandeo_real_antiguo=sigma_pandeo_real
                            end
                        end

                        # Si tiene una tension ingenieril de pandeo muy alta, la tension real puede hacerse casi infinito, así que limitamos su valor para poder efectuar los calculos 
                        if isinf(sigma_pandeo_real)==true; sigma_pandeo_real=sigma_real_limite_compresion; else ;end

                        epsilon_critico_pandeo_real= epsilon_02_real+(sigma_pandeo_real-sigma_02_real)/Young_02+epsilon_ult_p_asterisco*(sigma_pandeo_real-sigma_02_real)^m/(sigma_ultimo_real-sigma_02_real)^m
                        energia_critica_pandeo= sigma_pandeo_real*epsilon_critico_pandeo_real-(epsilon_02_real*sigma_pandeo_real+sigma_pandeo_real^2/(2*Young_02)-sigma_02_real*sigma_pandeo_real/Young_02+epsilon_ult_p_asterisco/(m+1)*(sigma_pandeo_real-sigma_02_real)^(m+1)/(sigma_ultimo_real-sigma_02_real)^m-epsilon_02_real*sigma_02_real+0.002/(n+1)*sigma_02_real+sigma_02_real^2/2*(1/Young+1/Young_02))
                    end             

                    sigma_ingenieril=sigma_ing_pandeo    
                    energias_especificas[elementos[i]]=sigma_pandeo_real*(epsilon_real-epsilon_critico_pandeo_real)+ energia_critica_pandeo
                    Youngs_tangentes_ing[elemento]=0.0001
                end
            end       

        end

    end

    return Youngs_tangentes_ing, energias_especificas

end #function curva_tension_deformacion 


function curva_tension_deformacion_tabla_Osgood(epsilons_ingenieril, nombre_invariante, tensiones_pandeo) # MODULO SECANTE --> obtener el modulo elastico equivalente (dañado) 

    ### Para calcular las energias especificas o los youngs secantes, hay que descomentarlos en lo de abajo
    Youngs_tangentes=Dict(); # energias_especificas=Dict(); Youngs_secantes=Dict()

    curva_tension_deformacion_tabla=readdlm(raw"curva_ingenieril_tension_deformacion_Osgood.txt") #leer el archivo para recoger datos del ensayo
    Young=Young_material()
    poisson=Poisson()


    deformacion=curva_tension_deformacion_tabla[1:end,1]
    tension=curva_tension_deformacion_tabla[1:end,2]
    energia=curva_tension_deformacion_tabla[1:end,3]
    Young_tangente=curva_tension_deformacion_tabla[1:end,4]



    ### BUSCAR EL VALOR DE EPSILON DE LA LISTA MAS CERCANO AL VALOR DE ENTRADA
    N_elementos=length(epsilons_ingenieril)
    elementos=collect(keys(epsilons_ingenieril))

    for i in 1:N_elementos

        if nombre_invariante=="inv_1"; epsilon_ingenieril=epsilons_ingenieril[elementos[i]]/(1-2*poisson)
        elseif nombre_invariante=="inv_2"; epsilon_ingenieril=epsilons_ingenieril[elementos[i]]/sqrt(2*poisson-poisson^2)
        elseif nombre_invariante=="inv_3"; epsilon_ingenieril=epsilons_ingenieril[elementos[i]]/cbrt(poisson^2)
        elseif nombre_invariante=="deform_equiv"; epsilon_ingenieril=epsilons_ingenieril[elementos[i]]
        elseif nombre_invariante=="deform_equiv_2"; epsilon_ingenieril=epsilons_ingenieril[elementos[i]]
        else println("no se ha introducido un criterio para entrar en la curva tension-deformacion")
        end
        
        # TRACCION
        if epsilon_ingenieril>=0 
            
            # DEFORMACION NULA
            if epsilon_ingenieril==0  # deformacion nula
                Youngs_tangentes[elementos[i]]=Young
                # energias_especificas[elementos[i]]=0.0               
                # Youngs_secantes[elementos[i]]=Young

            # PARTE ELASTICA + PLASTICA + FRACTURA (ABLANDAMIENTO)
            elseif epsilon_ingenieril < deformacion[end] # no se rompe el elemento
                j=1
                while deformacion[j]-epsilon_ingenieril<0; j=j+1; end
                sigma_ingenieril=(tension[j]+tension[j-1])/2
                Youngs_tangentes[elementos[i]]=(Young_tangente[j]+Young_tangente[j-1])/2
                # energias_especificas[elementos[i]]=(energia[j]+energia[j-1])/2
                # Youngs_secantes[elementos[i]]=sigma_ingenieril/epsilon_ingenieril

            # ROTURA
            else;   
                println("ATENCION: la barra ", elementos[i], " se ha roto" )
                Youngs_tangentes[elementos[i]]=0.001 #La barra a roto. Valor muy pequeño de modulo de Young pero que no sea cero para evitar fallos al invertir las matrices
                # energias_especificas[elementos[i]]=0.0
                # Youngs_secantes[elementos[i]]=0.001
            end

        # COMPRESION
        else;  
            epsilon_ingenieril=abs(epsilon_ingenieril)
            sigma_pandeo=tensiones_pandeo[elementos[i]]


            # PARTE ELASTICA + PLASTICA EN COMPRESION + PANDEO           
            j=1
            while (deformacion[j]-epsilon_ingenieril<0) && (j<length(tension)); j=j+1; end
            sigma_ingenieril=(tension[j]+tension[j-1])/2

            # COMPROBAR PANDEO
            if sigma_ingenieril < sigma_pandeo # No hay pandeo               
                Youngs_tangentes[elementos[i]]=(Young_tangente[j]+Young_tangente[j-1])/2
                # Youngs_secantes[elementos[i]]=abs(sigma_ingenieril/epsilon_ingenieril)
                # energias_especificas[elementos[i]]=(energia[j]+energia[j-1])/2
                # if j==length(tension) # para modelar los elementos que se pasan de la tension maxima pero no pandean 
                #     energias_especificas[elementos[i]]=energias_especificas[elementos[i]]+(epsilon_ingenieril-deformacion[j])*sigma_ingenieril
                # else
                # end


            else; # Hay pandeo
                println("pandeo")
                ### Determinar la deformacion a partir de la cual hay pandeo y la energia correspondiente         
                j=1
                while (tension[j]-sigma_pandeo < 0) && (j<length(tension));    j=j+1;  end               
                Youngs_tangentes[elementos[i]]=0.001
                # energia_critica_pandeo=(energia[j]+energia[j-1])/2  
                # deformacion_critica_pandeo=(deformacion[j]+deformacion[j-1])/2
                # sigma_ingenieril=sigma_pandeo
                # energias_especificas[elementos[i]]= energia_critica_pandeo + (epsilon_ingenieril-deformacion_critica_pandeo)*sigma_pandeo
                # Youngs_secantes[elementos[i]]=abs(sigma_pandeo/epsilon_ingenieril)
            end

        end
    end
    
    return Youngs_tangentes #, Youngs_secantes, energias_especificas

end #function curva_tension_deformacion



function curva_tension_deformacion_tabla_Osgood_1_elemento(epsilon_ingenieril, nombre_invariante, tension_pandeo) # MODULO SECANTE --> obtener el modulo elastico equivalente (dañado) 

    Youngs_tangentes=0; energias_especificas=0

    curva_tension_deformacion_tabla=readdlm(raw"curva_ingenieril_tension_deformacion_Osgood.txt") #leer el archivo para recoger datos del ensayo
    Young=Young_material()
    poisson=Poisson()


    deformacion=curva_tension_deformacion_tabla[1:end,1]
    tension=curva_tension_deformacion_tabla[1:end,2]
    energia=curva_tension_deformacion_tabla[1:end,3]
    Young_tangente=curva_tension_deformacion_tabla[1:end,4]



    ### BUSCAR EL VALOR DE EPSILON DE LA LISTA MAS CERCANO AL VALOR DE ENTRADA

    if nombre_invariante=="inv_1"; epsilon_ingenieril=epsilon_ingenieril/(1-2*poisson)
    elseif nombre_invariante=="inv_2"; epsilon_ingenieril=epsilon_ingenieril/sqrt(2*poisson-poisson^2)
    elseif nombre_invariante=="inv_3"; epsilon_ingenieril=epsilon_ingenieril/cbrt(poisson^2)
    elseif nombre_invariante=="deform_equiv"; epsilon_ingenieril=epsilon_ingenieril
    elseif nombre_invariante=="deform_equiv_2"; epsilon_ingenieril=epsilon_ingenieril
    else println("no se ha introducido un criterio para entrar en la curva tension-deformacion")
    end

    # TRACCION
    if epsilon_ingenieril>=0 
        
        # DEFORMACION NULA
        if epsilon_ingenieril==0  # deformacion nula
            energia_especifica=0.0
            Young_tangente=Young

        # PARTE ELASTICA + PLASTICA + FRACTURA (ABLANDAMIENTO)
        elseif epsilon_ingenieril < deformacion[end] # no se rompe el elemento
            j=1
            while deformacion[j]-epsilon_ingenieril<0; j=j+1; end
            #sigma_ingenieril=(tension[j]+tension[j-1])/2
            energia_especifica=(energia[j]+energia[j-1])/2
            Young_tangente=(Young_tangente[j]+Young_tangente[j-1])/2

        # ROTURA
        else;   
            Young_tangente=0.001 #La barra a roto. Valor muy pequeño de modulo de Young pero que no sea cero para evitar fallos al invertir las matrices
            energia_especifica=0.0
        end

    # COMPRESION
    else  
        epsilon_ingenieril=abs(epsilon_ingenieril)
        sigma_pandeo=tension_pandeo


        # PARTE ELASTICA + PLASTICA EN COMPRESION + PANDEO           
        j=1
        while (deformacion[j]-epsilon_ingenieril<0) && (j<length(tension)); j=j+1; end
        sigma_ingenieril=(tension[j]+tension[j-1])/2

        # COMPROBAR PANDEO
        if sigma_ingenieril < sigma_pandeo # No hay pandeo
            energia_especifica=(energia[j]+energia[j-1])/2
            Young_tangente=(Young_tangente[j]+Young_tangente[j-1])/2
            if j==length(tension) # para modelar los elementos que se pasan de la tension maxima pero no pandean 
                energia_especifica=energia_especifica+(epsilon_ingenieril-deformacion[j])*sigma_ingenieril
            else
            end


        else; # Hay pandeo
            ### Determinar la deformacion a partir de la cual hay pandeo y la energia correspondiente         
            j=1
            while (tension[j]-sigma_pandeo < 0) && (j<length(tension));    j=j+1;  end
            deformacion_critica_pandeo=(deformacion[j]+deformacion[j-1])/2
            energia_critica_pandeo=(energia[j]+energia[j-1])/2  

            sigma_ingenieril=sigma_pandeo
            energia_especifica= energia_critica_pandeo + (epsilon_ingenieril-deformacion_critica_pandeo)*sigma_pandeo
            Young_tangente=0.001
        end

    end


    return Young_tangente, energia_especifica

end #function curva_tension_deformacion


function generar_tabla_tension_deformacion_Osgood()
    # IMPORTANTE: TODO EL RATO DEFORMACIONES Y TENSIONES REALES (NO INGENIERILES)
    
    fichero_curva_real=open(raw"curva_real_tension_deformacion_Osgood.txt","w") #crear el archivo para recoger datos del ensayo
    fichero_curva_ingenieril=open(raw"curva_ingenieril_tension_deformacion_Osgood.txt","w") #crear el archivo para recoger datos del ensayo

    ### PARAMETROS DEL MATERIAL
    Young=194000  # MPa
    poisson=Poisson()
    sigma_02_real=445 # MPa
    sigma_ultimo_real=730 # MPa
    epsilon_ultimo_real=0.51
    n=4.85
    epsilon_fractura=0.53  # valor inventado pero tipico en aceros frente al epsilon_ultimo_real

    ### PARAMETROS DERIVADOS
    epsilon_02_real=sigma_02_real/Young+0.002
    e=sigma_02_real/Young
    m=1+3.5*sigma_02_real/sigma_ultimo_real
    Young_02=Young/(1+0.002*n/e) # modulo de young tangente (la derivada) en 0.2    
    epsilon_ult_ingenieril=exp(epsilon_ultimo_real)-1
    epsilon_fractura_ingenieril=exp(epsilon_fractura)-1



    integral_curva=0; epsilon_ingenieril_anterior=0; sigma_ingenieril_anterior=0; 
    epsilon_ingenieril=0; sigma_ingenieril=0; Young_tangente=Young
    inc_sigma=0.25 # incremento en tensiones reales para obtener la tabla
    inc_epsilon_daño=0.001

    # las curvas tension real-deformacion real son practicamente iguales en traccion que en compresion
    println(fichero_curva_ingenieril, epsilon_ingenieril," ",sigma_ingenieril," ",integral_curva," ",Young_tangente) # Para deformacion nula

    for sigma_real in inc_sigma:inc_sigma:sigma_02_real
        sigma_real_siguiente=sigma_real+inc_sigma
        epsilon_real=sigma_real/Young+0.002*(sigma_real/sigma_02_real)^n
        epsilon_real_siguiente=(sigma_real_siguiente)/Young+0.002*(sigma_real_siguiente/sigma_02_real)^n
        epsilon_ingenieril=exp(epsilon_real)-1
        epsilon_ingenieril_siguiente=exp(epsilon_real_siguiente)-1
        sigma_ingenieril=sigma_real/(1+epsilon_ingenieril)
        sigma_ingenieril_siguiente=(sigma_real_siguiente)/(1+epsilon_ingenieril_siguiente)


        # integral de la curva tension_deformacion para obtener la energía correspondiente a cada deformacion
        # la integral da el mismo resultado haciendola en real que en ingenieril
        integral_curva=integral_curva+(epsilon_ingenieril-epsilon_ingenieril_anterior)*(sigma_ingenieril+sigma_ingenieril_anterior)/2
        Young_tangente_atrasado=(sigma_ingenieril-sigma_ingenieril_anterior)/(epsilon_ingenieril-epsilon_ingenieril_anterior)
        Young_tangente_adelantado=(sigma_ingenieril_siguiente-sigma_ingenieril)/(epsilon_ingenieril_siguiente-epsilon_ingenieril)
        Young_tangente=(Young_tangente_adelantado+Young_tangente_atrasado)/2
        epsilon_ingenieril_anterior=epsilon_ingenieril
        sigma_ingenieril_anterior=sigma_ingenieril

    
        println(fichero_curva_real, epsilon_real," ",sigma_real," ",integral_curva)
        println(fichero_curva_ingenieril, epsilon_ingenieril," ",sigma_ingenieril," ",integral_curva," ",Young_tangente)
        
    end

    
    for sigma_real in sigma_02_real+inc_sigma:inc_sigma:sigma_ultimo_real-inc_sigma
        sigma_real_siguiente=sigma_real+inc_sigma
        epsilon_real=epsilon_02_real+(sigma_real-sigma_02_real)/Young_02+epsilon_ultimo_real*((sigma_real-sigma_02_real)/(sigma_ultimo_real-sigma_02_real))^m
        epsilon_real_siguiente=epsilon_02_real+(sigma_real_siguiente-sigma_02_real)/Young_02+epsilon_ultimo_real*((sigma_real_siguiente-sigma_02_real)/(sigma_ultimo_real-sigma_02_real))^m
        epsilon_ingenieril=exp(epsilon_real)-1
        epsilon_ingenieril_siguiente=exp(epsilon_real_siguiente)-1
        sigma_ingenieril=sigma_real/(1+epsilon_ingenieril)
        sigma_ingenieril_siguiente=(sigma_real_siguiente)/(1+epsilon_ingenieril_siguiente)

        if epsilon_real>=epsilon_ultimo_real # esto es necesario porque el sigma_ultimo_real y el epsilon_ultimo_real no se corresponden completamente con la expresion analitica
            break
        else
        end
        

        # integral de la curva tension_deformacion para obtener la energía correspondiente a cada deformacion
        integral_curva=integral_curva+(epsilon_ingenieril-epsilon_ingenieril_anterior)*(sigma_ingenieril+sigma_ingenieril_anterior)/2
        Young_tangente_atrasado=(sigma_ingenieril-sigma_ingenieril_anterior)/(epsilon_ingenieril-epsilon_ingenieril_anterior)
        Young_tangente_adelantado=(sigma_ingenieril_siguiente-sigma_ingenieril)/(epsilon_ingenieril_siguiente-epsilon_ingenieril)
        Young_tangente=(Young_tangente_adelantado+Young_tangente_atrasado)/2
        epsilon_ingenieril_anterior=epsilon_ingenieril
        sigma_ingenieril_anterior=sigma_ingenieril


        println(fichero_curva_real, epsilon_real," ",sigma_real," ",integral_curva)
        println(fichero_curva_ingenieril, epsilon_ingenieril," ",sigma_ingenieril," ",integral_curva," ",Young_tangente)

    end
    

    # PARTE DE FRACTURA (ABLANDAMIENTO) -- linea recta
    for epsilon_real in epsilon_ultimo_real+inc_epsilon_daño:inc_epsilon_daño:epsilon_fractura

        sigma_ingenieril_anterior=sigma_ingenieril
        epsilon_ingenieril_anterior=epsilon_ingenieril
        epsilon_ingenieril=exp(epsilon_real)-1
        sigma_ingenieril_ult=sigma_ultimo_real/(1+epsilon_ult_ingenieril)
        sigma_ingenieril=sigma_ingenieril_ult-(epsilon_ingenieril-epsilon_ult_ingenieril)/(epsilon_fractura_ingenieril-epsilon_ult_ingenieril)*sigma_ingenieril_ult
        sigma_real=sigma_ingenieril*(1+epsilon_ingenieril)

        sigma_real_siguiente=sigma_real+inc_sigma
        epsilon_real_siguiente=epsilon_real+inc_epsilon_daño
        epsilon_ingenieril_siguiente=exp(epsilon_real_siguiente)-1
        sigma_ingenieril_siguiente=sigma_ingenieril_ult-(epsilon_ingenieril_siguiente-epsilon_ult_ingenieril)/(epsilon_fractura_ingenieril-epsilon_ult_ingenieril)*sigma_ingenieril_ult
        

        # integral de la curva tension_deformacion para obtener la energía correspondiente a cada deformacion
        integral_curva=integral_curva+(epsilon_ingenieril-epsilon_ingenieril_anterior)*(sigma_ingenieril+sigma_ingenieril_anterior)/2
        Young_tangente_atrasado=(sigma_ingenieril-sigma_ingenieril_anterior)/(epsilon_ingenieril-epsilon_ingenieril_anterior)
        Young_tangente_adelantado=(sigma_ingenieril_siguiente-sigma_ingenieril)/(epsilon_ingenieril_siguiente-epsilon_ingenieril)
        Young_tangente=(Young_tangente_adelantado+Young_tangente_atrasado)/2
        epsilon_ingenieril_anterior=epsilon_ingenieril
        sigma_ingenieril_anterior=sigma_ingenieril


        println(fichero_curva_real, epsilon_real," ",sigma_real," ",integral_curva)
        println(fichero_curva_ingenieril, epsilon_ingenieril," ",sigma_ingenieril," ",integral_curva," ",Young_tangente)

    end

    close(fichero_curva_real)
    close(fichero_curva_ingenieril)

    return

end


function Poisson()
    poisson=0.3
    return poisson
end #function poisson()


function Young_material()
    young=194000
    return young
end # function Young


end # module


#$$$$
# using Plots
# ## PLOTEAR CURVA RAMBERG OSGOOD
# a=collect(0:0.001:0.85)
# plot(a,curva_tension_deformacion_tabla_Osgood.(a), title="secant Young's modulus",label="Ramberg Osgood",lw=3)
# plot!(a,curva_tension_deformacion.(a),label="approximated";marker=(:circle,1))
# xlabel!("engineering strain")
# ylabel!("Secant Young (MPa)")


# include(raw"curva_sigma_epsilon_material.jl"); using .curva_sigma_epsilon_material
# using Plots
# using DelimitedFiles
# generar_tabla_tension_deformacion_Osgood()
# curva_tension_deformacion_tabla=readdlm(raw"curva_ingenieril_tension_deformacion_Osgood.txt")
# deformacion=curva_tension_deformacion_tabla[1:end,1]
# tension=curva_tension_deformacion_tabla[1:end,2]
# plot(deformacion, tension)


#include(raw"curva_sigma_epsilon_material.jl"); using .curva_sigma_epsilon_material; using Plots; using DelimitedFiles;curva_tension_deformacion_tabla=readdlm(raw"curva_ingenieril_tension_deformacion_Osgood.txt");deformacion=curva_tension_deformacion_tabla[1:end,1];tension=curva_tension_deformacion_tabla[1:end,2]; plot(deformacion, tension)

