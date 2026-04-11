import numpy as np
from typing import List, Dict

class FractalBrancher:
    """
    Generador de ramificación fractal para Bio-Kidney AI 2026.
    Integra la Ley de Murray y geometría fractal 3D.
    """
    def __init__(self, murray_alpha: float = 3.0, min_radius_um: float = 6.0):
        self.alpha = murray_alpha
        self.min_radius = min_radius_um / 10000.0  # de micras a cm
        self.segments = []

    def calculate_child_radii(self, r_parent: float) -> float:
        # Ley de Murray: r_p^3 = r_h1^3 + r_h2^3 -> r_h = r_p / (2^(1/3))
        return r_parent / (2 ** (1 / self.alpha))

    def _rotate_vector(self, v, axis, angle):
        """Fórmula de rotación de Rodrigues para 3D."""
        v = np.array(v)
        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)
        return v * np.cos(angle) + np.cross(axis, v) * np.sin(angle) + \
               axis * np.dot(axis, v) * (1 - np.cos(angle))

    def grow_tree(self, start_node: np.ndarray, direction: np.ndarray, 
                  radius: float, depth: int, max_depth: int, length: float):
        """Genera la red frondosa de forma recursiva."""
        
        if depth >= max_depth or radius < self.min_radius:
            return

        child_radius = self.calculate_child_radii(radius)
        child_length = length * 0.8 
        
        dir_norm = direction / np.linalg.norm(direction)
        
        # Crear un vector ortogonal para definir el plano de bifurcación
        ortho = np.array([-dir_norm[1], dir_norm[0], 0])
        if np.allclose(ortho, 0): 
            ortho = np.array([0, -dir_norm[2], dir_norm[1]])
        
        # Ángulo de bifurcación (aprox 35 grados)
        angle = 0.6 
        
        dir_child1 = self._rotate_vector(dir_norm, ortho, angle)
        dir_child2 = self._rotate_vector(dir_norm, ortho, -angle)

        end1 = start_node + dir_child1 * child_length
        end2 = start_node + dir_child2 * child_length

        # Guardar en formato compatible con tu pipeline de datos
        self.segments.append({'p0': start_node.tolist(), 'p1': end1.tolist(), 'r': child_radius, 'lvl': depth})
        self.segments.append({'p0': start_node.tolist(), 'p1': end2.tolist(), 'r': child_radius, 'lvl': depth})

        # Recursión
        self.grow_tree(end1, dir_child1, child_radius, depth + 1, max_depth, child_length)
        self.grow_tree(end2, dir_child2, child_radius, depth + 1, max_depth, child_length)

if __name__ == "__main__":
    fb = FractalBrancher()
    # Punto de inicio (origen), Dirección (hacia arriba), Radio inicial (0.02cm), prof 0, max 5, largo inicial 0.1cm
    fb.grow_tree(np.array([0,0,0]), np.array([0,0,1]), 0.02, 0, 5, 0.1)
    print(f"Prueba completada: {len(fb.segments)} micro-segmentos generados.")
