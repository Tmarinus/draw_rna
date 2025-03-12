import matplotlib.pyplot as plt
from draw_rna.draw_utils import *
from nooverlap import Pusher

class mpl(object):
    def __init__(self, ax, fig=None, dpi=72):
        # create the file
        #self.dpi = dpi
        #plt.figure(figsize=(w/self.dpi,h/self.dpi), dpi=self.dpi)
        self.ax=ax
        self.fig=fig
        plt.sca(self.ax)

    def line(self, x1, y1, x2, y2, stroke, width=1, alpha=1, **kwargs):
        """"""
        # print 'Line (%s %s %s %s %s)' % (x1, y1, x2, y2, color)
        stroke = convert_color(stroke)

        return plt.plot([x1,x2], [y1,y2], linewidth=width, c=stroke, zorder=0, alpha=alpha, **kwargs)

    def circle(self, x, y, radius, fill, stroke, alpha=1, gid=None, **kwargs):
        fill = convert_color(fill)
        scatter = plt.scatter(x,y,s=1,zorder=0,alpha=alpha, **kwargs)
        circ = plt.Circle((x,y),radius=radius, alpha=alpha, color=fill, zorder=1, linewidth=0, **kwargs)
        circ.set_gid(gid)
        self.ax.add_artist(circ)
        return scatter, circ

    def text(self, x, y, size, fill, align, string, alpha=1, gid=None, check_position=False, **kwargs):
        fill = convert_color(fill)
        text = plt.text(x,y,string, fontsize=size, color=fill, zorder=2, horizontalalignment='center', verticalalignment='center', alpha=alpha, **kwargs)
        if check_position and self.fig:
            check_overlap(self.fig, self.ax)
        text.set_gid(gid)
        return text
    ## rotated
    #self.__out.write(' <text x="%d" y="%d" font-family="sans_serif" font-size="%d" fill="%s" text-anchor="%s" transform="rotate(180 %d,%d)">%s</text>' % (x-10,y+10,size,fill,align,x,y,str))

    def clean_up(self):
        self.ax.set_aspect('equal')
        plt.axis('off')
        l, r = plt.xlim()
        plt.xlim(l-20,r+20)
        b,t = plt.ylim()
        plt.ylim(b-20,t+20)




def check_overlap(fig, ax, max_x_shift = 0.3, max_y_shift = 0.3):
    """Pushes all text-elements in ax of figure to avoid overlap.

    The maximum push per iteration is max_x_shift and max_y_shift. These
    are factors of the width and height of the text-boxes.

    Parameters:
    -----------
    fig : matplotlib.figure.Figure
        The figure containing the axes.
    ax : matplotlib.axes.Axes
        The axes containing the text-elements.

    max_x_shift : float
    max_y_shift : float

    """
    # get all the positions from texts
    texts = ax.texts
    pusher = Pusher()
    old_positions = {}
    # idx = 0
    for text in texts:
        is_moveable = True
        if text.get_text().upper() in ['A', 'G', 'C', 'U', 'T', "5'", "3'"]:
            is_moveable = False
        r = fig.canvas.get_renderer()
        expand = (1.0, 1.0)
        ext = text.get_window_extent(r).expanded(*expand).transformed(ax.transData.inverted())
        position = text.get_position()
        # old_positions[idx] = position
        # idx += 1


        x0 = position[0]
        y0 = position[1]

        d_left = position[0] - ext.xmin
        d_right = ext.xmax - position[0]
        d_top = ext.ymax - position[1]
        d_bottom = position[1] - ext.ymin

        pusher.add_box(x0,y0,d_left,d_right,d_top,d_bottom, is_moveable)

    # push the boxes
    pusher.push_free(max_x_shift, max_y_shift)
    # import numpy as np
    # re-position the text objects
    idx = 0
    for text in texts:
        # if len(text.get_text()) < 4:
        # if text.get_text() in ['A', 'G', 'C', 'U', 'T', "5'", "3'"]:
        #     idx += 1
        #     continue
        position = pusher.get_position(idx)
        text.set_position(position)
        # if int(position[0]) != int(old_positions[idx][0]) and int(position[1]) != int(old_positions[idx][1]): #  round(position[0], 1) != round(old_positions[idx][0]) and round(position[1], 1) != round(old_positions[idx][1])
        # if np.allclose(position, old_positions[idx], atol=1e-2): #  round(position[0], 1) != round(old_positions[idx][0]) and round(position[1], 1) != round(old_positions[idx][1])
        #     return True
        idx += 1
    # return False
