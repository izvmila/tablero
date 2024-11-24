from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Blast import NCBIWWW,NCBIXML
import streamlit as st
from stmol import showmol
from stmol import * # pip install stmol==0.0.9 , pip install ipython_genutils
import matplotlib.pyplot as plt
from collections import Counter
from PIL import Image
import py3Dmol  # pip install py3Dmol==2.0.0.post2
import pandas as pd
import requests # # python -m pip install requests
import io  # Import the 'io' module for StringIO
from io import StringIO
from io import BytesIO
import streamlit.components.v1 as components


# Creamos una barra
st.sidebar.header("Proteínas Operaciones 🧬")
sidebar_render = st.sidebar.radio("Opciones : ",["Inicio", "Análisis de secuencia", "Parámetros de la estructura", "Secuencia de aminoácidos de proteínas", "Visualizador de proteínas"])

# Página principal
if sidebar_render == "Inicio":
    st.title('🧬 **Bioinformática: Análisis de Proteínas**')

    # Estilo de texto y colores
    st.markdown("""
    <style>
        .main-title {
            color: #4CAF50;  /* verde claro */
            font-size: 40px;
            font-weight: bold;
            text-align: center;
        }

        .text-block {
            color: #88dd9f;  /* verde más claro */
            font-size: 18px;
        }
        .team {
            font-style: bold;
            font-size: 16px;
            color: #7b9edd;  /* azul */
        }
    </style>
    """, unsafe_allow_html=True)

    # Título
    st.markdown('<div class="main-title">Bienvenido al Análisis de Proteínas</div>', unsafe_allow_html=True)

    # Descripción y subsecciones
    st.markdown("""
    <div class="text-block">
        Este tablero tiene el objetivo de facilitar el análisis y visualización de proteínas a partir de sus secuencias y estructuras. 
        Explora diferentes herramientas interactivas para estudiar sus propiedades y estructura. Las secciones disponibles son:
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    - **🔬 Análisis de secuencia**: Carga archivos FASTA y analiza las secuencias de proteínas. Extrae información relevante como la composición de aminoácidos y propiedades biofísicas.
    - **🧬 Parámetros de la estructura**: Calcula características estructurales, como el peso molecular, el punto isoeléctrico y la estabilidad de las proteínas, con un análisis detallado a nivel molecular.
    - **🔍 Secuencia de aminoácidos de proteínas**: Visualiza la secuencia y la proporción de átomos de diversas proteínas, con gráficos que permiten una mejor interpretación de sus características.
    - **🌐 Visualización 3D de proteínas**: Introduce un código PDB y explora la estructura tridimensional de proteínas en modelos interactivos. Personaliza la visualización y observa la estructura desde diferentes perspectivas.
    """, unsafe_allow_html=True)

    # Línea divisoria
    st.markdown("<hr style='border:1px solid #ccc;'/>", unsafe_allow_html=True)

    # Mensaje motivador
    st.markdown("""
    <div class="text-block">
        ¡Explora las herramientas del lado izquierdo y haz un análisis profundo de las proteínas que te interesen!
    </div>
    """, unsafe_allow_html=True)

    # Información del equipo
    st.markdown("<hr style='border:1px solid #ccc;'/>", unsafe_allow_html=True)
    st.markdown("<div class='team'>Equipo:</div>", unsafe_allow_html=True)
    st.markdown("""
    - **Camila García Rascón**
    - **Valeria Jara Salomón**
    """, unsafe_allow_html=True)

# Creamos Análisis de Secuencia
if sidebar_render == "Análisis de secuencia":
    st.title("🔬 Análisis de Secuencia")
    st.markdown("Sube tu secuencia de proteína en formato **FASTA** para analizarla ⬇️")

    # Función para leer y decodificar archivo FASTA
    def read_fasta_file(uploaded_file):
        fasta_file_content = uploaded_file.read().decode()
        return list(SeqIO.parse(StringIO(fasta_file_content), "fasta"))

    # Función para mostrar el contenido del archivo FASTA
    def display_fasta_file(sequence_contents):
        if len(sequence_contents) == 0:
            st.error("⚠️ Lo sentimos, no se encontraron secuencias en el archivo 😟")
        else:
            for i, record in enumerate(sequence_contents, start=1):
                with st.expander(f"🧬 Secuencia {i}: {record.id}"):
                    st.markdown(f"**🔖 ID:** `{record.id}`")
                    st.markdown(f"**🧾 Descripción:** {record.description}")
                    st.markdown(f"**🧪 Secuencia de aminoácidos:**")
                    st.code(str(record.seq), language="text")
                    st.markdown(f"**📏 Longitud de la secuencia:** `{len(record.seq)}`")
                st.divider()  # Línea divisoria entre secuencias

    # Subir archivo FASTA
    uploaded_file = st.file_uploader("📂 Sube tu archivo FASTA", type=["fasta"], help="Solo se admiten archivos con extensión .fasta")

    if uploaded_file is not None:
        with st.spinner("Procesando archivo... 🕒"):
            sequence_contents = read_fasta_file(uploaded_file)
            st.success("¡Archivo procesado exitosamente! 🎉", icon="✅")
        display_fasta_file(sequence_contents)
        st.info("Puedes copiar la secuencia para realizar otras operaciones.", icon="✨")


# Parámetros de la estructura de proteínas
if sidebar_render == "Parámetros de la estructura":
    st.title("🔬 Cálculos de parámetros de la estructura de proteínas")
    st.markdown("Introduce tu secuencia de proteína y ajusta los valores necesarios para analizar sus propiedades estructurales. 🌟")

    # Entrada de secuencia y nivel de pH
    sequence_input = st.text_area("✍️ Ingresa la secuencia de aminoácidos:")
    pH = st.number_input("🌡️ ¿Con qué nivel de pH deseas analizar tu proteína?", min_value=0.0, max_value=14.0, value=7.0, step=0.1)
    
    if st.button("⚡ ¡Calcular!"):
        if not sequence_input:
            st.error("Por favor, ingresa una secuencia para calcular sus propiedades.")
        else:
            # Análisis de la proteína
            sequence_reference = ProteinAnalysis(str(sequence_input))
            st.markdown("<h3 style='text-align: center; color: #4CAF50;'>✨ Propiedades calculadas de la proteína ✨</h3>", unsafe_allow_html=True)

            # Número de aminoácidos
            st.markdown("**1️⃣ Número de aminoácidos:**")
            st.info(f"🔢 **Valor:** `{sequence_reference.count_amino_acids()}`")

            # Peso molecular
            molecular_weight = round(sequence_reference.molecular_weight(), 2)
            st.markdown("**2️⃣ Peso molecular:**")
            st.info(f"⚖️ **Peso molecular:** `{molecular_weight} Da`")

            # Aromaticidad con barra de progreso
            aromaticity = round(sequence_reference.aromaticity(), 2)
            st.markdown("**3️⃣ Aromaticidad:**")
            st.info("Proporción de aminoácidos aromáticos en la proteína.")
            st.progress(min(int(aromaticity * 100), 100))  # Barra de progreso de 0 a 100%

            # Índice de inestabilidad con barra visual
            instability_index = round(sequence_reference.instability_index(), 2)
            stability = "La proteína es inestable" if instability_index >= 40 else "La proteína es estable"
            st.markdown("**4️⃣ Índice de inestabilidad:**")
            st.info(f"📉 **Valor:** `{instability_index}` - ⚖️ **Estabilidad:** {stability}")
            st.progress(min(int(instability_index), 100))  # Barra de progreso del índice de inestabilidad

            # Punto isoeléctrico
            isoelectric_point = round(sequence_reference.isoelectric_point(), 2)
            st.markdown("**5️⃣ Punto isoeléctrico (pI):**")
            st.info(f"⚡ **Valor:** `{isoelectric_point}`")

            # Fracciones de estructura secundaria con gráfico de barras
            st.markdown("**6️⃣ Estructura secundaria:**")
            secondary_structure = sequence_reference.secondary_structure_fraction()
            helix_fraction = round(secondary_structure[0] * 100, 2)
            turn_fraction = round(secondary_structure[1] * 100, 2)
            sheet_fraction = round(secondary_structure[2] * 100, 2)

            # Crear gráfico de barras
            structure_data = {
                "Tipo": ["Hélices", "Giros", "Láminas"],
                "Porcentaje (%)": [helix_fraction, turn_fraction, sheet_fraction],
            }
            st.bar_chart(structure_data)

            # Carga de la proteína a un pH específico con gráfico interactivo
            st.markdown("**7️⃣ Carga de la proteína a diferentes niveles de pH:**")
            ph_range = [x / 10 for x in range(0, 141)]  # Rango de pH de 0.0 a 14.0
            charge_values = [sequence_reference.charge_at_pH(ph) for ph in ph_range]
            
            # Crear un gráfico de línea interactivo
            charge_data = {"pH": ph_range, "Carga": charge_values}
            st.line_chart(charge_data)

            # Mostrar la carga específica al pH seleccionado
            charge = round(sequence_reference.charge_at_pH(pH), 3)
            st.info(f"⚡ **Carga a pH {pH}:** `{charge}`")


import requests  


# Definir las proteínas y sus secuencias
proteinas = {
    "Insulina": "MALWMRLLPLLALLALWGPDPAAAFGPGGPLALTLSSSINQEGASQSTSQPLNSRWQRPVEEQELLPCEDPQVP",  # Reemplaza con la secuencia real de insulina
    "Glucagon": "MKSIYFVAGLFVMLVQGSWQRSLQDTEEKSRSFSASQADPLSDPDQMNEDKRHSQGTFTSDYSKYLDSRRAQDFVQWLMNTKRNRNNIAKRHDEFERHAEGTFTSDVSSYLEGQAAKEFIAWLVKGRGRRDFPEEVAIVEELGRRHADGSFSDEMNTILDNLAARDFINWLIQTKITDRK",  # Reemplaza con la secuencia real de glucagón
    "Hemoglobina": "MTQTPYEVIGQERLYQLIDHFYSLVEQDNRINHLFPGDFAETARKQKQFLTQFLGGPDLYTQEHGHPMLRMRHLPFPIDDKAKEAWLENMHTAITHAQLPHGAGDYLYERLRLTANHMVNIEN",  # Reemplaza con la secuencia real de hemoglobina
    "Colageno": "MHPGLWLLLVTLCLTEELAAAGEKSYGKPCGGQDCSGSCQCFPEKGARGRPGPIGIQGPTGPQGFTGSTGLSGLKGERGFPGLLGPYGPKGDKGPMGVPGFLGINGIPGHPGQPGPRGPPGLDGCNGTQGAVGFPGPDGYPGLLGPPGLPGQKGSKGDPVLAPGSFKGMKGDPGLPGLDGITGPQGAPGFPGAVGPAGPPGLQGPPGPPGPLGPDGNMGLGFQGEKGVKGDVGLPGPAGPPPSTGELEFMGFPKGKKGSKGEPGPKGFPGISGPPGFPGLGTTGEKGEKGEKGIPGLPGPRGPMGSEGVQGPPGQQGKKGTLGFPGLNGFQGIEGQKGDIGLPGPDVFIDIDGAVISGNPGDPGVPGLPGLKGDEGIQGLRGPSGVPGLPALSGVPGALGPQGFPGLKGDQGNPGRTTIGAAGLPGRDGLPGPPGPPGPPSPEFETETLHNKESGFPGLRGEQGPKGNLGLKGIKGDSGFCACDGGVPNTGPPGEPGPPGPWGLIGLPGLKGARGDRGSGGAQGPAGAPGLVGPLGPSGPKGKKGEPILSTIQGMPGDRGDSGSQGFRGVIGEPGKDGVPGLPGLPGLPGDGGQGFPGEKGLPGLPGEKGHPGPPGLPGNGLPGLPGPRGLPGDKGKDGLPGQQGLPGSKGITLPCIIPGSYGPSGFPGTPGFPGPKGSRGLPGTPGQPGSSGSKGEPGSPGLVHLPELPGFPGPRGEKGLPGFPGLPGKDGLPGMIGSPGLPGSKGATGDIFGAENGAPGEQGLQGLTGHKGFLGDSGLPGLKGVHGKPGLLGPKGERGSPGTPGQVGQPGTPGSSGPYGIKGKSGLPGAPGFPGISGHPGKKGTRGKKGPPGSIVKKGLPGLKGLPGNPGLVGLKGSPGSPGVAGLPALSGPKGEKGSVGFVGFPGIPGLPGISGTRGLKGIPGSTGKMGPSGRAGTPGEKGDRGNPGPVGIPSPRRPMSNLWLKGDKGSQGSAGSNGFPGPRGDKGEAGRPGPPGLPGAPGLPGIIKGVSGKPGPPGFMGIRGLPGLKGSSGITGFPGMPGESGSQGIRGSPGLPGASGLPGLKGDNGQTVEISGSPGPKGQPGESGFKGTKGRDGLIGNIGFPGNKGEDGKVGVSGDVGLPGAPGFPGVAGMRGEPGLPGSSGHQGAIGPLGSPGLIGPKGFPGFPGLHGLNGLPGTKGTHGTPGPSITGVPGPAGLPGPKGEKGYPGIGIGAPGKPGLRGQKGDRGFPGLQGPAGLPGAPGISLPSLIAGQPGDPGRPGLDGERGRPGPAGPPGPPGPSSNQGDTGDPGFPGIPGFSGLPGELGLKGMRGEPGFMGTPGKVGPPGDPGFPGMKGKAGARGSSGLQGDPGQTPTAEAVQVPPGPLGLPGIDGIPGLTGDPGAQGPVGLQGSKGLPGIPGKDGPSGLPGPPGALGDPGLPGLQGPPGFEGAPGQQGPFGMPGMPGQSMRVGYTLVKHSQSEQVPPCPIGMSQLWVGYSLLFVEGQEKAHNQDLGFAGSCLPRFSTMPFIYCNINEVCHYARRNDKSYWLSTTAPIPMMPVSQTQIPQYISRCSVCEAPSQAIAVHSQDITIPQCPLGWRSLWIGYSFLMHTAAGAEGGGQSLVSPGSCLEDFRATPFIECSGARGTCHYFANKYSFWLTTVEERQQFGELPVSETLKAGQLHTRVSRCQVCMKSL",  # Reemplaza con la secuencia real de colágeno
}


# Función para calcular la proporción de átomos
def calcular_proporcion(proteina):
    secuencia = proteinas.get(proteina)
    if secuencia is None:
        return None
    aa_count = Counter(secuencia)
    return aa_count


# Interfaz en Secuencia de aminoácidos de proteínas
if sidebar_render == "Secuencia de aminoácidos de proteínas":
    st.title("🔍 Secuencia de Aminoácidos de Proteínas")
    st.markdown("Aquí puedes observar las secuencias de 4 proteínas diferentes y su proporción de átomos. 🌟")

    proteina_seleccionada = st.selectbox("Selecciona una proteína", list(proteinas.keys()))

    # Calcular la proporción de átomos
    aa_count = calcular_proporcion(proteina_seleccionada)

    if aa_count:
        st.write(f"Proporción de átomos en la proteína {proteina_seleccionada}:")
        st.write(aa_count)

        # Mostrar gráfico
        fig, ax = plt.subplots()
        ax.bar(aa_count.keys(), aa_count.values())
        st.pyplot(fig)
    


# Visualizador de proteínas
if sidebar_render == "Visualizador de proteínas":
    st.title("🧬 Visualizador de estructuras proteicas")
    st.markdown("Introduce el código PDB de tu proteína para observar su estructura en 3D. ¡Personaliza la visualización como prefieras! 🌟")

    # Enlace para buscar códigos PDB
    st.markdown(
        """
         **¿No conoces el código PDB?**  
        Puedes buscar proteínas en la base de datos [Protein Data Bank (PDB)](https://www.rcsb.org/) haciendo clic en este enlace. 🌐
        """
    )

    # Entrada de datos por parte del usuario
    ## Código PDB
    PDB_Code = st.text_input("🔑 Introduce el código PDB (4 caracteres):", value='9DIX', help="El código PDB representa una estructura específica de proteína almacenada en la base de datos Protein Data Bank.")
    ## Tipo de visualización
    DISPLAY_TYPE = st.selectbox("🎨 Elige el tipo de visualización:", ["línea", "dibujos animados (cartoon)", "cruz", "palitos (stick)", "esfera"], index=1)
    ## Residuo específico
    RESIDUE_FOCUS = st.text_input("🔍 Ingresa un residuo de aminoácidos para enfocar (Ej: ALA):", help="Puedes especificar un residuo (como ALA, LYS) para centrarte en él. Deja en blanco si no es necesario.")
    ## Animación de giro
    SPIN_ANIMATION = st.selectbox("🔄 ¿Quieres animación de giro en la estructura?", [True, False], index=0, help="Activa o desactiva la animación de giro para observar la estructura desde todas las perspectivas.")
    ## Color de fondo
    color_pal = st.color_picker("🌈 Elige el color de fondo:", "#f7f2f5", help="Selecciona el color de fondo para tu visualización 3D.")

    # Función para obtener información de la proteína desde PDB
    def fetch_pdb_metadata(pdb_code):
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_code.lower()}"
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()
        else:
            return None

    # Función para visualizar la proteína
    def PDB_VISUALISER(PDB_Code, color_pal, DISPLAY_TYPE, SPIN_ANIMATION, RESIDUE_FOCUS):
        # Convertir el tipo de visualización al formato esperado por py3Dmol
        display_mapping = {
            "línea": "line",
            "dibujos animados (cartoon)": "cartoon",
            "cruz": "cross",
            "palitos (stick)": "stick",
            "esfera": "sphere"
        }
        display_type_mapped = display_mapping[DISPLAY_TYPE]

        # Crear la visualización 3D
        xyzview = py3Dmol.view(query=f'pdb:{PDB_Code}')
        xyzview.setStyle({display_type_mapped: {'color': 'spectrum'}})
        xyzview.setBackgroundColor(color_pal)
        xyzview.zoomTo()
        xyzview.spin(SPIN_ANIMATION)

        if RESIDUE_FOCUS:
            # Centrar en el residuo especificado
            showmol(render_pdb_resn(viewer=xyzview, resn_lst=[RESIDUE_FOCUS]), height=500, width=800)
        else:
            # Mostrar la molécula completa
            showmol(xyzview, height=500, width=800)

    # Botón para generar la estructura
    if st.button("🔬 ¡Visualizar estructura!"):
        with st.spinner(f"⚙️ Obteniendo la estructura PDB para el código {PDB_Code}..."):
            try:
                # Obtener información de la proteína desde PDB
                metadata = fetch_pdb_metadata(PDB_Code)
                if metadata:
                    st.markdown("### ℹ️ Información sobre la proteína")
                    st.write(f"**Nombre de la proteína**: {metadata.get('struct', {}).get('title', 'No disponible')}")
                    st.write(f"**Fecha de liberación**: {metadata.get('rcsb_accession_info', {}).get('initial_release_date', 'No disponible')}")
                    st.write(f"**Método de determinación estructural**: {metadata.get('exptl', [{}])[0].get('method', 'No disponible')}")
                    st.write(f"**Resolución**: {metadata.get('rcsb_entry_info', {}).get('resolution_combined', ['No disponible'])[0]} Å")

                # Visualizar la estructura
                PDB_VISUALISER(PDB_Code, color_pal, DISPLAY_TYPE, SPIN_ANIMATION, RESIDUE_FOCUS)
                st.success(f"Estructura del código PDB `{PDB_Code}` generada con éxito. 🎉", icon="✅")
            except Exception as e:
                st.error(f"❌ Hubo un error al obtener la estructura o la información asociada. Verifica el código PDB e inténtalo nuevamente. \nError: {e}")
