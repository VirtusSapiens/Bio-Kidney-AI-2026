"""
Módulo: Analizador de Proyecto BioKidney-AI
------------------------------------------
Este módulo realiza un escaneo estructural y de metadatos del ecosistema BioKidney-AI.
Genera una base de conocimiento en formato Markdown (CONTEXTO_PROYECTO.md) que
facilita la orquestación del proyecto y el análisis arquitectónico.

Ingeniería de Software Aplicada:
- Gestión de rutas mediante pathlib para portabilidad.
- Patrón de diseño procedimental modular.
- Manejo robusto de excepciones y codificaciones.
- Extracción avanzada de docstrings e imports.

Autor: BioKidney-AI Team
Fecha: Marzo 2026
"""

import os
import datetime
import logging
from pathlib import Path
from typing import List, Dict, Optional

# Configuración de logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

class AnalizadorProyecto:
    """
    Clase encargada de inspeccionar la estructura del proyecto y generar
    documentación técnica automática.
    """
    
    def __init__(self, ruta_raiz: str = ".", archivo_salida: str = "CONTEXTO_PROYECTO.md"):
        self.ruta_raiz = Path(ruta_raiz).resolve()
        self.archivo_salida = archivo_salida
        self.extensiones_interes = ('.py', '.md', '.csv', '.txt')
        self.carpetas_ignoradas = {'.git', '__pycache__', 'env', 'node_modules'}

    def obtener_estructura_texto(self) -> str:
        """Genera una representación visual de la estructura de directorios."""
        lineas = ["\n## 1. Estructura de Directorios\n```text"]
        
        for path in sorted(self.ruta_raiz.rglob('*')):
            # Filtrar carpetas ignoradas
            if any(part in self.carpetas_ignoradas for part in path.parts):
                continue
                
            depth = len(path.relative_to(self.ruta_raiz).parts) - 1
            indent = ' ' * 4 * depth
            
            if path.is_dir():
                lineas.append(f"{indent}{path.name}/")
            else:
                lineas.append(f"{indent}    {path.name}")
                
        lineas.append("```\n")
        return "\n".join(lineas)

    def analizar_archivo_python(self, path: Path) -> List[str]:
        """Extrae metadatos específicos de archivos Python."""
        info = ["- **Tipo:** Script de Python"]
        try:
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                lineas = f.readlines()
                
            # Extraer Imports
            imports = [l.strip() for l in lineas if l.startswith(('import ', 'from ')) and not l.strip().startswith('#')]
            if imports:
                resumen_imports = ", ".join(sorted(list(set(imports)))[:8])
                info.append(f"- **Dependencias:** `{resumen_imports}{'...' if len(imports) > 8 else ''}`")
            
            # Extraer Docstring (Propósito)
            docstring = self._extraer_docstring(lineas)
            if docstring:
                info.append(f"- **Propósito:** {docstring[:250]}...")
                
        except Exception as e:
            info.append(f"- *Error en análisis:* {str(e)}")
        return info

    def _extraer_docstring(self, lineas: List[str]) -> str:
        """Lógica interna para extraer el primer bloque de comentarios/docstring."""
        doc = ""
        en_bloque = False
        for line in lineas[:40]:
            clean_line = line.strip()
            if '"""' in clean_line or "'''" in clean_line:
                if clean_line.count('"""') >= 2 or clean_line.count("'''") >= 2:
                    return clean_line.replace('"""', '').replace("'''", "").strip()
                en_bloque = not en_bloque
                continue
            if en_bloque:
                doc += clean_line + " "
        return doc.strip()

    def generar_inventario_binario(self) -> List[str]:
        """Cuenta y clasifica archivos binarios y activos multimedia."""
        lineas = ["## 3. Inventario de Activos Binarios\n"]
        tipos = {
            '.blend': 'Modelos Blender',
            '.png': 'Imágenes/Visualizaciones',
            '.pdf': 'Informes',
            '.mp4': 'Videos',
            '.obj': 'Mallas 3D',
            '.docx': 'Documentos Word'
        }
        
        for ext, desc in tipos.items():
            archivos = list(self.ruta_raiz.rglob(f"*{ext}"))
            if archivos:
                nombres = ", ".join(sorted([a.name for archivos in [archivos[:5]] for a in archivos]))
                lineas.append(f"#### {desc} ({ext})")
                lineas.append(f"- **Total:** {len(archivos)}")
                lineas.append(f"- **Muestras:** `{nombres}{'...' if len(archivos) > 5 else ''}`")
        
        return lineas

    def ejecutar(self):
        """Orquestador principal del análisis."""
        logging.info(f"Iniciando análisis profesional en: {self.ruta_raiz}")
        
        contenido = [
            f"# Informe de Contexto: BioKidney-AI",
            f"**Fecha de generación:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"**Estado del Sistema:** Pipeline 100% Completado (Hito Marzo 2026)",
            self.obtener_estructura_texto(),
            "## 2. Resumen de Componentes Clave\n"
        ]

        # Análisis de archivos de interés
        archivos_interes = sorted(
            [p for p in self.ruta_raiz.rglob('*') 
             if p.suffix in self.extensiones_interes and not any(part in self.carpetas_ignoradas for part in p.parts)],
            key=lambda x: str(x)
        )

        for path in archivos_interes:
            rel_path = path.relative_to(self.ruta_raiz)
            contenido.append(f"### Archivo: `{rel_path}`")
            
            if path.suffix == '.py':
                contenido.extend(self.analizar_archivo_python(path))
            elif path.suffix == '.md':
                contenido.append("- **Tipo:** Documentación Markdown")
            elif path.suffix == '.csv':
                try:
                    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                        header = f.readline().strip()
                        lineas = sum(1 for _ in f)
                        contenido.append("- **Tipo:** Dataset / Datos Estructurales")
                        contenido.append(f"- **Columnas:** `{header}`")
                        contenido.append(f"- **Registros:** {lineas}")
                except: pass
            
            contenido.append("\n---\n")

        # Activos binarios
        contenido.extend(self.generar_inventario_binario())

        # Guardar resultado
        try:
            with open(self.archivo_salida, 'w', encoding='utf-8') as f:
                f.write("\n".join(contenido))
            logging.info(f"Análisis completado exitosamente. Archivo: {self.archivo_salida}")
        except Exception as e:
            logging.error(f"Error al guardar el informe: {e}")

if __name__ == "__main__":
    analizador = AnalizadorProyecto()
    analizador.ejecutar()
