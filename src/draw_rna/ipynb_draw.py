import matplotlib.pyplot as plt
import draw_rna.draw as d
from draw_rna.draw_utils import seq2col

def draw_struct(seq, secstruct, c=None, ax=None, return_ax=False, show_plot=True, **kwargs):
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
        assert len(c) == len(secstruct), f"len(c): {len(c)} == len(secstruct): {len(secstruct)}"
        if isinstance(c[0], float):
            d.draw_rna(seq, secstruct, c, ext_color_file=True, ax=ax, show_plot=show_plot, **kwargs)
        else:
            d.draw_rna(seq, secstruct, c, ax=ax, show_plot=show_plot, **kwargs)

    else:
        d.draw_rna(seq, secstruct, seq2col(seq), ax=ax, show_plot=show_plot, **kwargs)
    # if scalar_map:
    #     plt.colorbar(scalar_map, label='Energies', ax=plt.gca())
    if return_ax:
        return plt.gca()
    if ax is None and show_plot:
        plt.show()
        print(show_plot)
    # plt.close('all')
    #return plt.gca()