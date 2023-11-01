#MODULO PARA GENERAR LOS FICHEROS PARA PLOTEAR EN PARAVIEW SOLO UN DESPLAZAMIENTO (VALE PARA ELEMENTOS BARRA, ELEMENTOS SOLIDOS DE DISTINTOS NUMEROS DE LADOS, etc.) 
module representar_paraview

using DelimitedFiles

export plot_paraview_estructura_metamaterial, plot_paraview_solido_tetraedros, plot_paraview_solido_hexaedros

function plot_paraview_estructura_metamaterial(coords,conectividad,desplazamientos,giros,fuerzas,momentos)
    N_nodos=length(coords[:,1])
    N_elementos=length(conectividad[:,1])

    fichero_plot=open(raw"representaciones_paraview_metamaterial.vtu","w") #crear el archivo
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
    
    n=1
    for i in 1:N_nodos
        println(fichero_plot, desplazamientos[n]," ",desplazamientos[n+1]," ",desplazamientos[n+2])
        n=n+3
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # GIROS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Giros(rad)" NumberOfComponents="3" format="ascii"> \n""")

    n=1
    for i in 1:N_nodos
        println(fichero_plot, giros[n]," ",giros[n+1]," ",giros[n+2])
        n=n+3
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # FUERZAS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Fuerzas(N)" NumberOfComponents="3" format="ascii"> \n""")
    
    n=1
    for i in 1:N_nodos
        println(fichero_plot, fuerzas[n]," ",fuerzas[n+1]," ",fuerzas[n+2])
        n=n+3
    end

    write(fichero_plot,"""              </DataArray> \n""")

    # MOMENTOS
    write(fichero_plot,"""              <DataArray type="Float64" Name="Momentos(N·mm)" NumberOfComponents="3" format="ascii"> \n""")

    n=1
    for i in 1:N_nodos
        println(fichero_plot, momentos[n]," ",momentos[n+1]," ",momentos[n+2])
        n=n+3
    end

    write(fichero_plot,"""              </DataArray> \n""")

    write(fichero_plot,"""          </PointData> \n""")
    write(fichero_plot,"""       </Piece> \n""")
    write(fichero_plot,"""   </UnstructuredGrid> \n""")
    write(fichero_plot,"""</VTKFile> \n""")

    close(fichero_plot)
end # function plot_paraview_estructura_metamaterial


function plot_paraview_solido_tetraedros(coords, conectividad, u, f)
    N_nodos=length(coords[:,1])
    N_elementos=length(conectividad[:,1])

    fichero_plot=open(raw"representaciones_paraview_solido_tetraedros.vtu","w") #crear el archivo
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


function plot_paraview_solido_hexaedros(coords, conectividad, u, f)
    N_nodos=length(coords[:,1])
    N_elementos=length(conectividad[:,1])

    fichero_plot=open(raw"representaciones_paraview_solido_hexaedros.vtu","w") #crear el archivo
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


end # module