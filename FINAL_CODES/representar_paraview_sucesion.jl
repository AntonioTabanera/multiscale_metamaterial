#MODULO PARA GENERAR LOS FICHEROS PARA PLOTEAR EN PARAVIEW CUANDO HAY UNA SUCESIÓN DE DESPLAZAMIENTOS HASTA EL VALOR FINAL
module representar_paraview_sucesion

using DelimitedFiles

export plot_paraview_estructura_metamaterial_sucesion,plot_paraview_solido_sucesion_hexaedros, plot_paraview_estructura_metamaterial_multiescala_sucesion_hexaedros, plot_paraview_estructura_metamaterial_multiescala_sucesion_tetraedros, plot_paraview_solido_sucesion_tetraedros

function plot_paraview_estructura_metamaterial_sucesion(coords,conectividad,desplazamientos,giros,fuerzas,momentos, iter)
    N_nodos=length(coords[:,1])
    N_elementos=length(conectividad[:,1])

    fichero_plot=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\representaciones_paraview\\3D\\metamaterial\\metamaterial_$(lpad(iter,2,"0")).vtu","w") #crear el archivo
    write(fichero_plot, """<?xml version="1.0"?> \n""")

    ### GEOMETRIA
    write(fichero_plot, """<VTKFile type="UnstructuredGrid"> \n""")
    write(fichero_plot,"""   <UnstructuredGrid> \n""")
    println(fichero_plot,"""      <Piece NumberOfPoints=" """, N_nodos ,""" "  NumberOfCells=" """, N_elementos ,""" ">""")
    write(fichero_plot,"""          <Points> \n""")
    write(fichero_plot,"""             <DataArray type="Float64" NumberOfComponents="3" format="ascii"> \n""")
    
    for i in 1:N_nodos
        println(fichero_plot, coords[i,1]," ",coords[i,2]," ",coords[i,3])  # Las coordenadas de los puntos
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""          </Points> \n""")
    write(fichero_plot,"""          <Cells> \n""")
    write(fichero_plot,"""             <DataArray type="Int32" Name="connectivity" format="ascii"> \n""")    

    for i in 1:N_elementos
        println(fichero_plot, conectividad[i,1]-1," ",conectividad[i,2]-1)  # La conectividad de los puntos
    end

    write(fichero_plot,"""              </DataArray> \n""")
    write(fichero_plot,"""              <DataArray type="Int32" Name="offsets" format="ascii"> \n""")

    for i in 1:N_elementos
        println(fichero_plot, 2*i)            # salto de nodos que componen los elementos (obligatorio para paraview)
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""              <DataArray type="UInt8" Name="types" format="ascii"> \n""")

    for i in 1:N_elementos
        println(fichero_plot, 3)            # tipo de elementos (segun paraview) --> lineas=3
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""          </Cells> \n""")  
    
    ### VALORES DE VARIABLES EN LOS NODOS
    write(fichero_plot,"""          <PointData> \n""")

    # DESPLAZAMIENTOS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Displacements(mm)" NumberOfComponents="3" format="ascii"> \n""")
    
    for nodo in coords[:,4]
        println(fichero_plot, desplazamientos[nodo][1]," ",desplazamientos[nodo][2]," ",desplazamientos[nodo][3])
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # GIROS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Giros(rad)" NumberOfComponents="3" format="ascii"> \n""")

    for nodo in coords[:,4]
        println(fichero_plot, giros[nodo][1]," ",giros[nodo][2]," ",giros[nodo][3])
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # FUERZAS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Fuerzas(N)" NumberOfComponents="3" format="ascii"> \n""")
    
    for nodo in coords[:,4]
        println(fichero_plot, fuerzas[nodo][1]," ",fuerzas[nodo][2]," ",fuerzas[nodo][3])
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # MOMENTOS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Momentos(N·mm)" NumberOfComponents="3" format="ascii"> \n""")

    for nodo in coords[:,4]
        println(fichero_plot, momentos[nodo][1]," ",momentos[nodo][2]," ",momentos[nodo][3])
    end

    write(fichero_plot,"""              </DataArray> \n""")

    write(fichero_plot,"""          </PointData> \n""")
    write(fichero_plot,"""       </Piece> \n""")
    write(fichero_plot,"""   </UnstructuredGrid> \n""")
    write(fichero_plot,"""</VTKFile> \n""")

    close(fichero_plot)
end # function plot_paraview_estructura_metamaterial


function plot_paraview_solido_sucesion_tetraedros(coords, conectividad, u, f, iter)
    N_nodos=length(coords[:,1])
    N_elementos=length(conectividad[:,1])

    fichero_plot=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\representaciones_paraview\\3D\\solido_tetraedros\\solido_tetraedros_$(lpad(iter,2,"0")).vtu","w") #crear el archivo
    write(fichero_plot, """<?xml version="1.0"?> \n""")

    ### GEOMETRIA
    write(fichero_plot, """<VTKFile type="UnstructuredGrid"> \n""")
    write(fichero_plot,"""   <UnstructuredGrid> \n""")
    println(fichero_plot,"""      <Piece NumberOfPoints=" """, N_nodos ,""" "  NumberOfCells=" """, N_elementos ,""" ">""")
    write(fichero_plot,"""          <Points> \n""")
    write(fichero_plot,"""             <DataArray type="Float64" NumberOfComponents="3" format="ascii"> \n""")
    
    for i in 1:N_nodos
        println(fichero_plot, coords[i,1]," ",coords[i,2]," ",coords[i,3])  # Las coordenadas de los puntos
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""          </Points> \n""")
    write(fichero_plot,"""          <Cells> \n""")
    write(fichero_plot,"""             <DataArray type="Int32" Name="connectivity" format="ascii"> \n""")    

    for i in 1:N_elementos
        println(fichero_plot, conectividad[i,1]-1," ",conectividad[i,2]-1, " ",conectividad[i,3]-1, " ",conectividad[i,4]-1)  # La conectividad de los puntos
    end

    write(fichero_plot,"""              </DataArray> \n""")
    write(fichero_plot,"""              <DataArray type="Int32" Name="offsets" format="ascii"> \n""")

    for i in 1:N_elementos
        println(fichero_plot, 4*i)            # salto de nodos que componen los elementos (obligatorio para paraview)
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""              <DataArray type="UInt8" Name="types" format="ascii"> \n""")

    for i in 1:N_elementos
        println(fichero_plot, 10)            # tipo de elementos (segun paraview) --> tertaedros=10
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""          </Cells> \n""")  
    
    ### VALORES DE VARIABLES EN LOS NODOS
    write(fichero_plot,"""          <PointData> \n""")

    # DESPLAZAMIENTOS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Displacements(mm)" NumberOfComponents="3" format="ascii"> \n""")
    
    n=1
    for i in 1:N_nodos
        println(fichero_plot, u[n]," ",u[n+1]," ",u[n+2])
        n=n+3
    end

    write(fichero_plot,"""              </DataArray> \n""")


    # FUERZAS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Fuerzas(N)" NumberOfComponents="3" format="ascii"> \n""")
    
    n=1
    for i in 1:N_nodos
        println(fichero_plot, f[n]," ",f[n+1]," ",f[n+2])
        n=n+3
    end

    write(fichero_plot,"""              </DataArray> \n""")


    write(fichero_plot,"""          </PointData> \n""")
    write(fichero_plot,"""       </Piece> \n""")
    write(fichero_plot,"""   </UnstructuredGrid> \n""")
    write(fichero_plot,"""</VTKFile> \n""")

    close(fichero_plot)
end # function plot_paraview_estructura_solida


function plot_paraview_solido_sucesion_hexaedros(coords, conectividad, u, f, iter)
    N_nodos=length(coords[:,1])
    N_elementos=length(conectividad[:,1])

    fichero_plot=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\representaciones_paraview\\3D\\solido_hexaedros\\solido_hexaedros_$(lpad(iter,2,"0")).vtu","w") #crear el archivo
    write(fichero_plot, """<?xml version="1.0"?> \n""")

    ### GEOMETRIA
    write(fichero_plot, """<VTKFile type="UnstructuredGrid"> \n""")
    write(fichero_plot,"""   <UnstructuredGrid> \n""")
    println(fichero_plot,"""      <Piece NumberOfPoints=" """, N_nodos ,""" "  NumberOfCells=" """, N_elementos ,""" ">""")
    write(fichero_plot,"""          <Points> \n""")
    write(fichero_plot,"""             <DataArray type="Float64" NumberOfComponents="3" format="ascii"> \n""")
    
    for i in 1:N_nodos
        println(fichero_plot, coords[i,1]," ",coords[i,2]," ",coords[i,3])  # Las coordenadas de los puntos
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""          </Points> \n""")
    write(fichero_plot,"""          <Cells> \n""")
    write(fichero_plot,"""             <DataArray type="Int32" Name="connectivity" format="ascii"> \n""")    

    for i in 1:N_elementos
        println(fichero_plot, conectividad[i,1]-1," ",conectividad[i,3]-1, " ",conectividad[i,4]-1, " ",conectividad[i,2]-1, " ",conectividad[i,5]-1," ",conectividad[i,7]-1, " ",conectividad[i,8]-1, " ",conectividad[i,6]-1)  # La conectividad de los puntos
    end

    write(fichero_plot,"""              </DataArray> \n""")
    write(fichero_plot,"""              <DataArray type="Int32" Name="offsets" format="ascii"> \n""")

    for i in 1:N_elementos
        println(fichero_plot, 8*i)            # salto de nodos que componen los elementos (obligatorio para paraview)
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""              <DataArray type="UInt8" Name="types" format="ascii"> \n""")

    for i in 1:N_elementos
        println(fichero_plot, 12)            # tipo de elementos (segun paraview) --> hexaedros=12
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""          </Cells> \n""")  
    
    ### VALORES DE VARIABLES EN LOS NODOS
    write(fichero_plot,"""          <PointData> \n""")

    # DESPLAZAMIENTOS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Displacements(mm)" NumberOfComponents="3" format="ascii"> \n""")
    
    n=1
    for i in 1:N_nodos
        println(fichero_plot, u[n]," ",u[n+1]," ",u[n+2])
        n=n+3
    end

    write(fichero_plot,"""              </DataArray> \n""")


    # FUERZAS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Fuerzas(N)" NumberOfComponents="3" format="ascii"> \n""")
    
    n=1
    for i in 1:N_nodos
        println(fichero_plot, f[n]," ",f[n+1]," ",f[n+2])
        n=n+3
    end

    write(fichero_plot,"""              </DataArray> \n""")


    write(fichero_plot,"""          </PointData> \n""")
    write(fichero_plot,"""       </Piece> \n""")
    write(fichero_plot,"""   </UnstructuredGrid> \n""")
    write(fichero_plot,"""</VTKFile> \n""")

    close(fichero_plot)
end # function plot_paraview_estructura_solida


function plot_paraview_estructura_metamaterial_multiescala_sucesion_tetraedros(coords_zona_metam, conectividad_metam_tetraed, desplazamientos, giros, fuerzas, momentos, elemento, iter)
    N_nodos=length(coords_zona_metam[:,1])
    N_elementos=length(conectividad_metam_tetraed[:,1])


    ### GENERAR Y ESCRIBIR EL FICHERO
    fichero_plot=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\representaciones_paraview\\3D\\multiescala_tetraedros\\metamaterial_elemento$(lpad(elemento,2,"0"))_$(lpad(iter,2,"0")).vtu","w") #crear el archivo
    write(fichero_plot, """<?xml version="1.0"?> \n""")

    ### GEOMETRIA
    write(fichero_plot, """<VTKFile type="UnstructuredGrid"> \n""")
    write(fichero_plot,"""   <UnstructuredGrid> \n""")
    println(fichero_plot,"""      <Piece NumberOfPoints=" """, N_nodos ,""" "  NumberOfCells=" """, N_elementos ,""" ">""")
    write(fichero_plot,"""          <Points> \n""")
    write(fichero_plot,"""             <DataArray type="Float64" NumberOfComponents="3" format="ascii"> \n""")
    
    for i in 1:N_nodos
        println(fichero_plot, coords_zona_metam[i,1]," ",coords_zona_metam[i,2]," ",coords_zona_metam[i,3])  # Las coordenadas de los puntos
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""          </Points> \n""")
    write(fichero_plot,"""          <Cells> \n""")
    write(fichero_plot,"""             <DataArray type="Int32" Name="connectivity" format="ascii"> \n""")    

    for i in 1:N_elementos
        nodo_1=findall(in(conectividad_metam_tetraed[i,1]),coords_zona_metam[:,4])[1]
        nodo_2=findall(in(conectividad_metam_tetraed[i,2]),coords_zona_metam[:,4])[1]
        println(fichero_plot, nodo_1-1," ",nodo_2-1)  # La conectividad de los puntos
    end

    write(fichero_plot,"""              </DataArray> \n""")
    write(fichero_plot,"""              <DataArray type="Int32" Name="offsets" format="ascii"> \n""")

    for i in 1:N_elementos
        println(fichero_plot, 2*i)            # salto de nodos que componen los elementos (obligatorio para paraview)
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""              <DataArray type="UInt8" Name="types" format="ascii"> \n""")

    for i in 1:N_elementos
        println(fichero_plot, 3)            # tipo de elementos (segun paraview) --> lineas=3
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""          </Cells> \n""")  
    
    ### VALORES DE VARIABLES EN LOS NODOS
    write(fichero_plot,"""          <PointData> \n""")

    # DESPLAZAMIENTOS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Displacements(mm)" NumberOfComponents="3" format="ascii"> \n""")
    

    for nodo in coords_zona_metam[:,4]
        println(fichero_plot, desplazamientos[nodo][1]," ",desplazamientos[nodo][2]," ",desplazamientos[nodo][3])
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # GIROS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Giros(rad)" NumberOfComponents="3" format="ascii"> \n""")

    
    for nodo in coords_zona_metam[:,4]
        println(fichero_plot, giros[nodo][1]," ",giros[nodo][2]," ",giros[nodo][3])
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # FUERZAS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Fuerzas(N)" NumberOfComponents="3" format="ascii"> \n""")
    
    n=1
    for nodo in coords_zona_metam[:,4]
        println(fichero_plot, fuerzas[nodo][1]," ",fuerzas[nodo][2]," ",fuerzas[nodo][3])
        n=n+3
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # MOMENTOS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Momentos(N·mm)" NumberOfComponents="3" format="ascii"> \n""")

   
    for nodo in coords_zona_metam[:,4]
        println(fichero_plot, momentos[nodo][1]," ",momentos[nodo][2]," ",momentos[nodo][3])
    end

    write(fichero_plot,"""              </DataArray> \n""")

    write(fichero_plot,"""          </PointData> \n""")
    write(fichero_plot,"""       </Piece> \n""")
    write(fichero_plot,"""   </UnstructuredGrid> \n""")
    write(fichero_plot,"""</VTKFile> \n""")

    close(fichero_plot)
end # function plot_paraview_estructura_tetraedro_metamaterial


function plot_paraview_estructura_metamaterial_multiescala_sucesion_hexaedros(coords_zona_metam, conectividad_metam_tetraed, desplazamientos, giros, fuerzas, momentos, elemento, iter)
    N_nodos=length(coords_zona_metam[:,1])
    N_elementos=length(conectividad_metam_tetraed[:,1])


    ### GENERAR Y ESCRIBIR EL FICHERO
    fichero_plot=open("C:\\__UNIVERSIDAD__\\AEROESPACIAL\\4º GIA\\TFG\\representaciones_paraview\\3D\\multiescala_hexaedros\\metamaterial_elemento$(lpad(elemento,2,"0"))_$(lpad(iter,2,"0")).vtu","w") #crear el archivo
    write(fichero_plot, """<?xml version="1.0"?> \n""")

    ### GEOMETRIA
    write(fichero_plot, """<VTKFile type="UnstructuredGrid"> \n""")
    write(fichero_plot,"""   <UnstructuredGrid> \n""")
    println(fichero_plot,"""      <Piece NumberOfPoints=" """, N_nodos ,""" "  NumberOfCells=" """, N_elementos ,""" ">""")
    write(fichero_plot,"""          <Points> \n""")
    write(fichero_plot,"""             <DataArray type="Float64" NumberOfComponents="3" format="ascii"> \n""")
    
    for i in 1:N_nodos
        println(fichero_plot, coords_zona_metam[i,1]," ",coords_zona_metam[i,2]," ",coords_zona_metam[i,3])  # Las coordenadas de los puntos
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""          </Points> \n""")
    write(fichero_plot,"""          <Cells> \n""")
    write(fichero_plot,"""             <DataArray type="Int32" Name="connectivity" format="ascii"> \n""")    

    for i in 1:N_elementos
        nodo_1=findall(in(conectividad_metam_tetraed[i,1]),coords_zona_metam[:,4])[1]
        nodo_2=findall(in(conectividad_metam_tetraed[i,2]),coords_zona_metam[:,4])[1]
        println(fichero_plot, nodo_1-1," ",nodo_2-1)  # La conectividad de los puntos
    end

    write(fichero_plot,"""              </DataArray> \n""")
    write(fichero_plot,"""              <DataArray type="Int32" Name="offsets" format="ascii"> \n""")

    for i in 1:N_elementos
        println(fichero_plot, 2*i)            # salto de nodos que componen los elementos (obligatorio para paraview)
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""              <DataArray type="UInt8" Name="types" format="ascii"> \n""")

    for i in 1:N_elementos
        println(fichero_plot, 3)            # tipo de elementos (segun paraview) --> lineas=3
    end

    write(fichero_plot,"""             </DataArray> \n""")
    write(fichero_plot,"""          </Cells> \n""")  
    
    ### VALORES DE VARIABLES EN LOS NODOS
    write(fichero_plot,"""          <PointData> \n""")

    # DESPLAZAMIENTOS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Displacements(mm)" NumberOfComponents="3" format="ascii"> \n""")
    

    for nodo in coords_zona_metam[:,4]
        println(fichero_plot, desplazamientos[nodo][1]," ",desplazamientos[nodo][2]," ",desplazamientos[nodo][3])
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # GIROS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Giros(rad)" NumberOfComponents="3" format="ascii"> \n""")

    
    for nodo in coords_zona_metam[:,4]
        println(fichero_plot, giros[nodo][1]," ",giros[nodo][2]," ",giros[nodo][3])
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # FUERZAS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Fuerzas(N)" NumberOfComponents="3" format="ascii"> \n""")
    
    n=1
    for nodo in coords_zona_metam[:,4]
        println(fichero_plot, fuerzas[nodo][1]," ",fuerzas[nodo][2]," ",fuerzas[nodo][3])
        n=n+3
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # MOMENTOS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Momentos(N·mm)" NumberOfComponents="3" format="ascii"> \n""")

   
    for nodo in coords_zona_metam[:,4]
        println(fichero_plot, momentos[nodo][1]," ",momentos[nodo][2]," ",momentos[nodo][3])
    end

    write(fichero_plot,"""              </DataArray> \n""")

    write(fichero_plot,"""          </PointData> \n""")
    write(fichero_plot,"""       </Piece> \n""")
    write(fichero_plot,"""   </UnstructuredGrid> \n""")
    write(fichero_plot,"""</VTKFile> \n""")

    close(fichero_plot)
end # function plot_paraview_estructura_tetraedro_metamaterial

end # module