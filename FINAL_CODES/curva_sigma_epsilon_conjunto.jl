#MODULO PARA ELMODULO DE YOUNG Y DE POISON DEL SOLIDO CON PROPIEDADES HOMOGENEAS (SOLO LA PARTE ELASTICA DE LA CURVA)
module curva_sigma_epsilon_conjunto

export curva_tension_deformacion_conjunto, Poisson_conjunto

# ESTOS VALORES DEL YOUNG Y DEL POISSON SE DEBEN AJUSTAR CON VALORES EXPERIMENTALES DE UNA PROBETA DE METAMATERIAL
# EL YOUNG ES SOLO PARA LA ZONA ELASTICA (SIN DAÑO)
function curva_tension_deformacion_conjunto(N_divisiones) # obtener el modulo elastico equivalente (dañado)
    if N_divisiones==10
        Young_medio_metam=433584.9902779429#50000#494056.37# 10          50000#761220.9#462378.4539939174#813522.9#760441#929498#813544.6#152417.4596807553#491108#
    elseif N_divisiones==11
        Young_medio_metam=488414.55# 11
    elseif N_divisiones==12
        Young_medio_metam=483732.515# 12
    elseif N_divisiones==13
        Young_medio_metam=479784.388# 13
    elseif N_divisiones==14
        Young_medio_metam=476410.307# 14
    elseif N_divisiones==15
        Young_medio_metam=473493.633# 15
    elseif N_divisiones==16
        Young_medio_metam=470839.739# 16
    elseif N_divisiones==17
        Young_medio_metam=468704.99# 17
    elseif N_divisiones==18
        Young_medio_metam=466715.369# 18
    elseif N_divisiones==19
        Young_medio_metam=464937.99# 19
    else
        println("numero de divisiones no introducido")
        Young_medio_metam=464937.99
    end

    return Young_medio_metam
end #function curva_tension_deformacion 


function Poisson_conjunto()
    poisson=0.13#389# 78926665# 0.14186#0.2699738898898168 #0.3373#0.187275892 0.19832093#
    return poisson
end #function poisson()


end # module