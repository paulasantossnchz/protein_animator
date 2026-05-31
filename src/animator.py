# animator.py
import os
import subprocess

def generar_animacion_pymol(lista_rutas_pdbs, ruta_salida_gif):
    exito_animacion = False
    
    if len(lista_rutas_pdbs) > 1:
        ruta_script_temporal = "script_pymol_temp.py"
        
        # 1. Construimos el código que PyMOL ejecutará en su propio proceso aislado
        lineas_script = [
            "import os",
            "from pymol import cmd",
            "cmd.reinitialize()"
        ]
        
        for indice in range(len(lista_rutas_pdbs)):
            ruta_segura = lista_rutas_pdbs[indice].replace("\\", "/")
            lineas_script.append("cmd.load('" + ruta_segura + "', 'estr_" + str(indice) + "')")
            
            if indice > 0:
                lineas_script.append("cmd.align('estr_" + str(indice) + "', 'estr_0')")
                
            lineas_script.append("cmd.create('trayectoria_completa', 'estr_" + str(indice) + "', 1, " + str(indice + 1) + ")")
            
        lineas_script.append("cmd.morph('animacion_final', 'trayectoria_completa')")
        lineas_script.append("fotogramas = cmd.count_states('animacion_final')")
        lineas_script.append("cmd.mset('1 x' + str(fotogramas))")
        
        # Configuramos la calidad y guardamos
        ruta_gif_segura = ruta_salida_gif.replace("\\", "/")
        lineas_script.append("cmd.save('" + ruta_gif_segura + "', 'animacion_final')")
        lineas_script.append("cmd.quit()")
        
        # Convertimos la lista de líneas en un único string de texto
        contenido_script = "\n".join(lineas_script)
        
        # 2. Escribimos el script temporal en el disco (aplicando la regla de apertura segura)
        archivo_script = open(ruta_script_temporal, 'w')
        try:
            archivo_script.write(contenido_script)
            escritura_correcta = True
        except Exception:
            escritura_correcta = False
        finally:
            archivo_script.close()
            
# 3. Ejecutamos PyMOL como un proceso independiente si se escribió el script
        if escritura_correcta:
            try:
                # Lanzamos PyMOL y capturamos la salida de texto (capture_output=True, text=True)
                proceso = subprocess.run(["pymol", "-c", "-q", ruta_script_temporal], capture_output=True, text=True)
                
                # 4. Comprobamos si el subproceso logró generar el GIF físico
                if os.path.exists(ruta_salida_gif):
                    exito_animacion = True
                else:
                    # Imprimimos el error real de PyMOL en tu terminal de Linux
                    print("\n=== LOG DE ERROR DE PYMOL ===")
                    print("SALIDA ESTÁNDAR:", proceso.stdout)
                    print("ERRORES:", proceso.stderr)
                    print("=============================\n")
                    exito_animacion = False
            except Exception as e:
                print("Error de Python al intentar llamar a PyMOL:", str(e))
                exito_animacion = False
            finally:
                # 5. Limpiamos el rastro eliminando el script temporal
                if os.path.exists(ruta_script_temporal):
                    os.remove(ruta_script_temporal)