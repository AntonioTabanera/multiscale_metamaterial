using DelimitedFiles
using Plots

Nx=20

X=10*Nx; vol=X^3
N_inc=50; inc_energia_total=zeros(N_inc); inc_desplaz=0.05#1.5e0/N_inc
 inc_deformaciones=inc_desplaz/X

inc_energias_elementos=readdlm(raw"C:\__UNIVERSIDAD__\AEROESPACIAL\4º GIA\TFG\COMPARACIONES\comparaciones_metamaterial_inc_energ_total_TRACC_2.txt")

for i in 1:N_inc
    inc_energia_total[i]=inc_energias_elementos[i]#sum(inc_energias_elementos[1+64*(i-1):64*i])
end


Youngs_tangentes=inc_energia_total/(1/2*inc_deformaciones^2*vol)


fichero_curva_real=open(raw"young_tangente_vs_deformacion_solido_homogeneo.txt","w")

epsilon_acumulado=0; sigma_acumulado=0
for young_tangente in Youngs_tangentes

    global epsilon_acumulado, inc_deformaciones, sigma_acumulado

    inc_sigma=young_tangente*inc_deformaciones

    epsilon_acumulado=epsilon_acumulado+inc_deformaciones
    sigma_acumulado=sigma_acumulado+inc_sigma
    println(fichero_curva_real, young_tangente," ", epsilon_acumulado," ", sigma_acumulado)
    
end

close(fichero_curva_real)



#### REPRESENTAR CURVA SIGMA-EPSILON DEL SOLIDO
datos=readdlm(raw"young_tangente_vs_deformacion_solido_homogeneo.txt")
plot(datos[:,2], datos[:,3], label="σ-ϵ curve", legend=:bottomright)
xlabel!("engineering strain")
ylabel!("engineering stress [MPa]")
titulo=string("STRESS-STRAIN CURVE OF THE SOLID")

display(title!(titulo))