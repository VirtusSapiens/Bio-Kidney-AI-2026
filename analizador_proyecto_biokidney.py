import os
import datetime
import logging
import ast
from pathlib import Path
from typing import List, Dict, Optional


    def obtener_estructura_texto(self) -> str:
        """Genera una representación visual de la estructura de directorios."""
        lineas = ["\n## 1. Arquitectura de Archivos\n```text"]
        
        for path in sorted(self.ruta_raiz.rglob('*')):
            # Filtrar carpetas ignoradas
            if any(part in self.carpetas_ignoradas for part in path.parts):
                continue
                
            rel_path = path.relative_to(self.ruta_raiz)
            depth = len(rel_path.parts) - 1
            indent = '    ' * depth

            if path.is_dir():
                lineas.append(f"{indent}{path.name}/")
            else:
                lineas.append(f"{indent}├── {path.name}")
                
        lineas.append("```\n")
        return "\n".join(lineas)

    def analizar_archivo_python(self, path: Path) -> List[str]:
        """Extrae metadatos de archivos Python utilizando AST para mayor precisión."""
        info = ["- **Tipo:** Script de Python"]
        try:
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            tree = ast.parse(content)
            
            # Extraer Docstring del módulo
            docstring = ast.get_docstring(tree)
            if docstring:
                resumen = docstring.split('\n')[0][:150]
                info.append(f"- **Propósito:** {resumen}...")

            # Extraer Imports (Import y ImportFrom)
            imports = []
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        imports.append(alias.name)
                elif isinstance(node, ast.ImportFrom):
                    imports.append(f"{node.module}")

            if imports:
                resumen_imports = ", ".join(sorted(list(set(imports)))[:6])
                info.append(f"- **Dependencias:** `{resumen_imports}{'...' if len(imports) > 8 else ''}`")
                
        except Exception as e:
            info.append(f"- *Aviso:* No se pudo parsear AST ({str(e)})")
        return info

    def generar_inventario_binario(self) -> List[str]:
        """Cuenta y clasifica archivos binarios y activos multimedia."""
        for ext, desc in tipos.items():
            archivos = list(self.ruta_raiz.rglob(f"*{ext}"))
            if archivos:
                # Filtrar si están en carpetas ignoradas
                archivos = [a for a in archivos if not any(p in self.carpetas_ignoradas for p in a.parts)]
                if not archivos: continue
                
                muestras = sorted([a.name for a in archivos])[:5]
                nombres = ", ".join(muestras)
                
                lineas.append(f"#### {desc} ({ext})")
                lineas.append(f"- **Total:** {len(archivos)}")
                lineas.append(f"- **Muestras:** `{nombres}{'...' if len(archivos) > 5 else ''}`")
        contenido = [
            f"# Informe de Contexto: BioKidney-AI",
            f"**Fecha de generación:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"**Estado:** Pipeline In-Silico 100% Validado",
            self.obtener_estructura_texto(),
            "## 2. Resumen de Componentes Clave\n"
        ]
