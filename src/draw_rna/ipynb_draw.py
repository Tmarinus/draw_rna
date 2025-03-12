import matplotlib.pyplot as plt
import draw_rna.draw as d
from draw_rna.draw_utils import seq2col

def draw_struct(seq, secstruct, c=None, scalar_map=None, line=False, large_mode=False, store_plt_svg=False,
 cmap='viridis', rotation=0, vmin=None, vmax=None, alpha=None, ax=None, numbering=None, custom_text=None, svg_mode=False, filename="secstruct",
                return_ax=False, show_plot=True,custom_line=False,
                custom_line_col='#6d7075', custom_line_width=2, custom_line_alpha=0.8, custom_line_style='dashed'):
    '''
    Draw sequence with secondary structure.
    Inputs:
    c (string or array-like).  If string, characters must correspond to colors.
     If array-like obj, used as mapping for colormap (provided in cmap), or a string.
    line (bool): draw secstruct as line.
    large_mode: draw outer loop as straight line.
    rotation: rotate molecule (in degrees).
    '''
    if seq == None:
        c = ['black' for x in secstruct]
    if c is not None:
        assert len(c) == len(secstruct)
        if isinstance(c[0], float):
            d.draw_rna(seq, secstruct, c, color_map=scalar_map, line=line, ext_color_file=True, cmap_name = cmap, vmin=vmin, vmax=vmax, store_plt_svg=store_plt_svg,
            rotation=rotation, large_mode = large_mode, alpha=alpha, ax=ax, numbering=numbering, custom_text=custom_text, svg_mode=svg_mode, filename=filename, show_plot=show_plot, custom_line=custom_line,
            custom_line_col=custom_line_col, custom_line_width=custom_line_width, custom_line_alpha=custom_line_alpha, custom_line_style=custom_line_style)
        else:
            d.draw_rna(seq, secstruct, c, color_map=scalar_map,  line=line, cmap_name=cmap, large_mode=large_mode, vmin=vmin, vmax=vmax, store_plt_svg=store_plt_svg,
            rotation=rotation, alpha=alpha, ax=ax, numbering=numbering, custom_text=custom_text, svg_mode=svg_mode, filename=filename, show_plot=show_plot, custom_line=custom_line,
            custom_line_col=custom_line_col, custom_line_width=custom_line_width, custom_line_alpha=custom_line_alpha, custom_line_style=custom_line_style)

    else:
        d.draw_rna(seq, secstruct, seq2col(seq),  color_map=scalar_map, line=line, cmap_name = cmap, vmin=vmin, vmax=vmax, store_plt_svg=store_plt_svg,
        large_mode = large_mode, rotation=rotation, alpha=alpha, ax=ax, numbering=numbering, custom_text=custom_text, svg_mode=svg_mode, filename=filename, show_plot=show_plot, custom_line=custom_line,
        custom_line_col=custom_line_col, custom_line_width=custom_line_width, custom_line_alpha=custom_line_alpha, custom_line_style=custom_line_style)
    # if scalar_map:
    #     plt.colorbar(scalar_map, label='Energies', ax=plt.gca())
    if return_ax:
        return plt.gca()
    if ax is None and show_plot:
        plt.show()
    # plt.close('all')
    #return plt.gca()