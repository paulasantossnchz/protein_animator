# app.py
import os
import streamlit as st
import py3Dmol
from stmol import showmol

# Importamos nuestros propios módulos
from src import pdb_fetcher
from src import animator
    

def leer_archivo_texto(ruta_archivo):
    contenido = ""
    archivo = open(ruta_archivo, 'r')
    try:
        contenido = archivo.read()
    finally:
        archivo.close()
    return contenido

def leer_archivo_binario(ruta_archivo):
    contenido = b""
    archivo = open(ruta_archivo, 'rb')
    try:
        contenido = archivo.read()
    finally:
        archivo.close()
    return contenido

def ejecutar_aplicacion():
    # Variable de estado para mantener el único punto de salida
    estado_ejecucion = True
    
    # 1. Título y descripción
    st.set_page_config(page_title="Animador de Proteínas", layout="centered")
    st.title("🧬 Animador de Proteínas en 3D")
    st.write("Introduce una lista de códigos PDB de la misma proteína para visualizar sus cambios conformacionales.")
    
    # 2. Input del usuario
    entrada_pdbs = st.text_input("Códigos PDB (separados por comas):", "1CRN, 1YYF")
    
    # 3. Botón de ejecución
    if st.button("Generar Animación"):
        if entrada_pdbs.strip() != "":
            # Limpiamos la entrada usando listas por comprensión (sin break/continue)
            lista_ids = []
            fragmentos = entrada_pdbs.split(",")
            for fragmento in fragmentos:
                id_limpio = fragmento.strip().upper()
                if id_limpio != "":
                    lista_ids.append(id_limpio)
            
            if len(lista_ids) > 1:
                st.info("📥 Descargando estructuras desde RCSB PDB...")
                carpeta_datos = os.path.join("data", "input_pdbs")
                ruta_salida_gif = os.path.join("data", "output_animations", "animacion.gif")
                
                # Aseguramos que la carpeta de salida exista
                if not os.path.exists(os.path.dirname(ruta_salida_gif)):
                    os.makedirs(os.path.dirname(ruta_salida_gif))
                
                # Llamada al módulo fetcher
                rutas_descargadas = pdb_fetcher.descargar_pdbs(lista_ids, carpeta_datos)
                
                if len(rutas_descargadas) > 1:
                    st.info("⚙️ Generando interpolación y renderizando con PyMOL (esto puede tardar)...")
                    
                    # Llamada al módulo animator
                    exito = animator.generar_animacion_pymol(rutas_descargadas, ruta_salida_gif)
                    
                    if exito:
                        st.success("¡Animación generada con éxito!")
                        
                        # 4. Mostrar el GIF usando nuestra función de lectura segura
                        datos_gif = leer_archivo_binario(ruta_salida_gif)
                        st.image(datos_gif, caption="Interpolación de los estados conformacionales")
                        
                        st.divider()
                        
                        # 5. Mostrar el primer PDB en 3D interactivo
                        st.subheader("Estructura de referencia (Interactiva)")
                        st.write(f"Mostrando conformación base: **{lista_ids[0]}**")
                        
                        contenido_pdb = leer_archivo_texto(rutas_descargadas[0])
                        
                        visor = py3Dmol.view(width=700, height=500)
                        visor.addModel(contenido_pdb, 'pdb')
                        
                        # Coloreado por estructura secundaria usando ssJmol
                        visor.setStyle({'cartoon': {'colorscheme': 'ssJmol'}})
                        visor.zoomTo()
                        
                        showmol(visor, height=500, width=700)
                        
                    else:
                        st.error("Error al generar la animación. Revisa que los PDBs sean de la misma proteína y tengan residuos equivalentes.")
                else:
                    st.error("No se pudieron descargar suficientes archivos. Comprueba los códigos PDB.")
            else:
                st.warning("Se necesitan al menos dos códigos PDB para crear una animación.")
        else:
            st.warning("Por favor, introduce al menos dos códigos PDB separados por comas.")

    # Único return al final de la función
    return estado_ejecucion

# Lanzamos la aplicación si es el script principal
if __name__ == "__main__":
    ejecutar_aplicacion()
