from dataclasses import dataclass, field
from typing import Optional, List, Tuple, Any, Union
import matplotlib.pyplot as plt

@dataclass
class VisualizationSettings:
    # Basic visual settings
    black_white: bool = False
    show_energies: bool = False
    energies_individual: bool = False
    store_plt_svg: bool = True
    line: bool = False
    large_mode: bool = False
    movie_mode: bool = False
    svg_mode: bool = False
    show_plot: bool = True
    
    # Colors and color mapping
    colors: Optional[List[Any]] = None
    color_map: Optional[str] = None
    cmap_name: str = 'viridis'
    vmin: Optional[float] = None
    vmax: Optional[float] = None
    alpha: Optional[List[float]] = None
    
    # Custom line settings
    custom_line: bool = False
    custom_line_col: str = '#6d7075'
    custom_line_width: float = 2
    custom_line_alpha: float = 0.8
    custom_line_style: str = 'dashed'

@dataclass
class StructureSettings:
    # Structure parameters
    predict: bool = False
    constraints: Optional[List[Any]] = None
    segment_lengths: bool = True
    numbering: Optional[List[Any]] = None
    add_fake: List[Any] = field(default_factory=list)
    
    # RNA rendering parameters
    NODE_R: float = 10
    PRIMARY_SPACE: float = 20
    PAIR_SPACE: float = 20
    external_multiplier: float = 1
    external_offset: float = 0
    rotation: float = 0

@dataclass
class TextSettings:
    # Text elements
    custom_text: List[Any] = field(default_factory=list)
    title: Optional[Tuple[str, float, float]] = None
    energy_title: bool = True
    RENDER_IN_LETTERS: bool = True
    text_size: float = 15

@dataclass
class FileSettings:
    # File related settings
    file_path: Optional[str] = None
    file_name: Optional[str] = None
    seq: Optional[str] = None
    db: Optional[str] = None

@dataclass
class PlotSettings:
    # Main categories as nested dataclasses
    visual: VisualizationSettings = field(default_factory=VisualizationSettings)
    structure: StructureSettings = field(default_factory=StructureSettings)
    text: TextSettings = field(default_factory=TextSettings)
    file: FileSettings = field(default_factory=FileSettings)
    
    # Matplotlib specific settings
    ax: Optional[plt.Axes] = None
    fig: Optional[plt.Figure] = None
    
    # Cell padding and size settings
    CELL_PADDING: int = 40
    cell_size_x: Optional[float] = None
    cell_size_y: Optional[float] = None

    def __post_init__(self):
        # Adjust settings based on mode
        if self.visual.large_mode or self.visual.movie_mode:
            self.CELL_PADDING = 100
            self.structure.external_multiplier = 10000
            self.structure.external_offset = 3.14159 + 3.14159 * self.structure.rotation/180
        else:
            self.CELL_PADDING = 40
            self.structure.external_multiplier = 1
            self.structure.external_offset = 0 + 3.14159 * self.structure.rotation/180

# Usage example:
settings = PlotSettings(
    visual=VisualizationSettings(
        black_white=True,
        show_energies=True,
        custom_line=True
    ),
    structure=StructureSettings(
        predict=False,
        segment_lengths=True,
        rotation=45
    ),
    text=TextSettings(
        title=("RNA Structure", 0.5, 0.95),
        energy_title=True
    ),
    file=FileSettings(
        file_name="output"
    )
)