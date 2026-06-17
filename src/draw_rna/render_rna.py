import draw_rna.svg as svg
import re, random, math
import numpy as np
from nooverlap import Pusher
import sys

class RNATreeNode:

    def __init__(self):
        self.children_ = []
        self.is_pair_ = False
        self.index_a_ = -1
        self.index_b_ = -1
        self.x_ = 0
        self.y_ = 0
        self.go_x_ = 0
        self.go_y_ = 0

def get_pairmap_from_secstruct(secstruct):
    """
    generates dictionary containing pair mappings

    args:
    secstruct contains secondary structure string

    returns:
    dictionary with pair mappings
    """
    pair_stack = []
    end_stack = []
    pairs_array = []
    i_range = list(range(0,len(secstruct)))

    # initialize all values to -1, meaning no pair
    for ii in i_range:
        pairs_array.append(-1)

    # assign pairs based on secstruct
    for ii in i_range:
        if(secstruct[ii] == "("):
            pair_stack.append(ii)
        elif(secstruct[ii] == ")"):
            if not pair_stack:
                end_stack.append(ii)
            else:
                index = pair_stack.pop()
                pairs_array[index] = ii
                pairs_array[ii] = index
    if len(pair_stack) == len(end_stack):
        n = len(pair_stack)
        for ii in range(n):
            pairs_array[pair_stack[ii]] = end_stack[-ii]
            pairs_array[end_stack[-ii]] = pair_stack[ii]
    else:
         print("ERROR: pairing incorrect %s" % secstruct)

    return pairs_array


def add_nodes_recursive(bi_pairs, rootnode, start_index, end_index):

    if(start_index > end_index) :
        print("Error occured while drawing RNA %d %d" % (start_index, end_index))
        sys.exit(0)

    if(bi_pairs[start_index] == end_index) :

        newnode = RNATreeNode()
        newnode.is_pair_ = True
        newnode.index_a_ = start_index
        newnode.index_b_ = end_index

        add_nodes_recursive(bi_pairs, newnode, start_index+1, end_index-1)

    else :
        newnode = RNATreeNode()
        jj = start_index
        while jj <= end_index:
            if(bi_pairs[jj] > jj) :
                add_nodes_recursive(bi_pairs,newnode, jj, bi_pairs[jj])
                jj = bi_pairs[jj] + 1
            else :
                newsubnode = RNATreeNode()
                newsubnode.is_pair_ = False
                newsubnode.index_a_ = jj
                newnode.children_.append(newsubnode)
                jj += 1

    rootnode.children_.append(newnode)

def setup_coords_recursive(rootnode, parentnode, start_x, start_y, go_x, go_y, NODE_R, PRIMARY_SPACE, PAIR_SPACE, external_multiplier, external_offset):

    cross_x = -go_y
    cross_y = go_x
    # print(rootnode.index_a_, rootnode.index_b_)
    children_width = len(rootnode.children_) * NODE_R * 2

    #print('children_width', children_width)

    rootnode.go_x_ = go_x
    rootnode.go_y_ = go_y

    if(len(rootnode.children_) == 1):
        rootnode.x_ = start_x
        rootnode.y_ = start_y

        if(rootnode.children_[0].is_pair_):
            setup_coords_recursive(rootnode.children_[0], rootnode, start_x + go_x * PRIMARY_SPACE, start_y + go_y * PRIMARY_SPACE, go_x, go_y, NODE_R, PRIMARY_SPACE, PAIR_SPACE, external_multiplier, external_offset)
        elif(rootnode.children_[0].is_pair_ == False and rootnode.children_[0].index_a_ < 0):
            setup_coords_recursive(rootnode.children_[0], rootnode, start_x, start_y, go_x, go_y, NODE_R, PRIMARY_SPACE, PAIR_SPACE, external_multiplier, external_offset)
        else:
            setup_coords_recursive(rootnode.children_[0], rootnode, start_x + go_x * PRIMARY_SPACE, start_y + go_y * PRIMARY_SPACE, go_x, go_y, NODE_R, PRIMARY_SPACE, PAIR_SPACE, external_multiplier, external_offset)

    elif(len(rootnode.children_) > 1) :

        npairs = 0
        for ii in range(0, len(rootnode.children_)):
            if(rootnode.children_[ii].is_pair_) :
                npairs+=1

        circle_length = (len(rootnode.children_) + 1) * PRIMARY_SPACE + (npairs + 1) * PAIR_SPACE
        circle_radius = circle_length / (2 * math.pi)

        length_walker = PAIR_SPACE / 2.0

        if (parentnode == None) :
            rootnode.x_ = go_x * circle_radius
            rootnode.y_ = go_y * circle_radius
            circle_radius *= external_multiplier
        else :
            rootnode.x_ = parentnode.x_ + go_x * circle_radius
            rootnode.y_ = parentnode.y_ + go_y * circle_radius

        for ii in range(0,len(rootnode.children_)):

            if (parentnode ==None):

                length_walker += PRIMARY_SPACE

                if(rootnode.children_[ii].is_pair_) :
                    length_walker += PAIR_SPACE / 2.0

                rad_angle = length_walker/circle_length * 2 * math.pi / external_multiplier - math.pi / 2.0 + external_offset
                child_x = rootnode.x_ + math.cos(rad_angle) * cross_x * circle_radius + math.sin(rad_angle) * go_x * circle_radius
                child_y = rootnode.y_ + math.cos(rad_angle) * cross_y * circle_radius + math.sin(rad_angle) * go_y * circle_radius

                child_go_x = child_x - rootnode.x_
                child_go_y = child_y - rootnode.y_
                child_go_len = math.sqrt(child_go_x * child_go_x + child_go_y * child_go_y)

                setup_coords_recursive(rootnode.children_[ii], rootnode, child_x, child_y, child_go_x / child_go_len, child_go_y / child_go_len, NODE_R, PRIMARY_SPACE, PAIR_SPACE, external_multiplier, external_offset)

                if(rootnode.children_[ii].is_pair_) :
                    length_walker += PAIR_SPACE / 2.0

            else:

                length_walker += PRIMARY_SPACE

                if(rootnode.children_[ii].is_pair_) :
                    length_walker += PAIR_SPACE / 2.0

                rad_angle = length_walker/circle_length * 2 * math.pi - math.pi / 2.0
                child_x = rootnode.x_ + math.cos(rad_angle) * cross_x * circle_radius + math.sin(rad_angle) * go_x * circle_radius
                child_y = rootnode.y_ + math.cos(rad_angle) * cross_y * circle_radius + math.sin(rad_angle) * go_y * circle_radius

                child_go_x = child_x - rootnode.x_
                child_go_y = child_y - rootnode.y_
                child_go_len = math.sqrt(child_go_x * child_go_x + child_go_y * child_go_y)

                setup_coords_recursive(rootnode.children_[ii], rootnode, child_x, child_y, child_go_x / child_go_len, child_go_y / child_go_len, NODE_R, PRIMARY_SPACE, PAIR_SPACE, external_multiplier, external_offset)

                if(rootnode.children_[ii].is_pair_) :
                    length_walker += PAIR_SPACE / 2.0

    else :
        rootnode.x_ = start_x
        rootnode.y_ = start_y

def get_coords_recursive(rootnode, xarray, yarray, PRIMARY_SPACE, PAIR_SPACE):
    if(rootnode.is_pair_) :
        cross_x = -rootnode.go_y_
        cross_y = rootnode.go_x_

        xarray[rootnode.index_a_] = rootnode.x_ + cross_x * PAIR_SPACE/2.0
        xarray[rootnode.index_b_] = rootnode.x_ - cross_x * PAIR_SPACE/2.0

        yarray[rootnode.index_a_] = rootnode.y_ + cross_y * PAIR_SPACE/2.0
        yarray[rootnode.index_b_] = rootnode.y_ - cross_y * PAIR_SPACE/2.0
    elif(rootnode.index_a_ >= 0) :
        xarray[rootnode.index_a_] = rootnode.x_
        yarray[rootnode.index_a_] = rootnode.y_

    for ii in range(0, len(rootnode.children_)):
        get_coords_recursive(rootnode.children_[ii], xarray, yarray, PRIMARY_SPACE, PAIR_SPACE)



class RNARenderer:

    def __init__(self):
        self.root_ = None
        self.xarray_ = None
        self.yarray_ = None
        self.size_ = None

    def setup_tree(self, secstruct, NODE_R,PRIMARY_SPACE, PAIR_SPACE, external_multiplier, external_offset):

        dangling_start = 0
        dangling_end = 0
        bi_pairs = get_pairmap_from_secstruct(secstruct)

        self.NODE_R = NODE_R
        self.root_ = None

        for ii in range(0,len(bi_pairs)):
            if bi_pairs[ii] < 0:
                dangling_start+=1
            else:
                break

        for ii in (len(bi_pairs)-1,-1, -1) :
            if(bi_pairs[ii] < 0):
                dangling_end+=1
            else:
                break

        self.root_ = RNATreeNode()

        #for jj in range(0,len(bi_pairs)):
        jj = 0
        while (jj < len(bi_pairs)):
            if (bi_pairs[jj] > jj) :
                add_nodes_recursive(bi_pairs,self.root_, jj, bi_pairs[jj])
                jj = bi_pairs[jj] + 1
            else:
                newsubnode = RNATreeNode()
                newsubnode.is_pair_ = False
                newsubnode.index_a_ = jj
                self.root_.children_.append(newsubnode)
                jj += 1
        xarray = []
        yarray = []

        for ii in range(0,len(secstruct)):
            xarray.append(0.0)
            yarray.append(0.0)

        self.setup_coords(NODE_R,PRIMARY_SPACE,PAIR_SPACE, external_multiplier, external_offset)
        self.get_coords(xarray,yarray,PRIMARY_SPACE,PAIR_SPACE)

        min_x = xarray[0] - NODE_R
        min_y = yarray[0] - NODE_R
        max_x = xarray[0] + NODE_R
        max_y = xarray[0] + NODE_R

        for x in xarray:
            if x - NODE_R < min_x:
                min_x = x - NODE_R
            if x + NODE_R > max_x:
                max_x = x + NODE_R

        for y in yarray:
            if y - NODE_R < min_y:
                min_y = y - NODE_R
            if y + NODE_R > max_y:
                max_y = y + NODE_R

        for ii in range(0,len(xarray)):
            xarray[ii] -= min_x
            yarray[ii] -= min_y

        self.size_ = [max_x - min_x, max_y - min_y]
        self.xarray_ = xarray
        self.yarray_ = yarray

    def get_size(self):
        return self.size_

    def draw(self, svgobj, offset_x, offset_y, colors, pairs, sequence, render_in_letter, external_offset, line=False, svg_mode=True, alpha=None,
             numbering=None, custom_text=None, custom_line=False, show_primes=True,
             custom_line_col='#6d7075', custom_line_width=2, custom_line_alpha=0.8, custom_line_style='dashed', **kwargs):
        if alpha is None:
            alpha = np.ones(len(self.xarray_))

        if self.xarray_ != None:

            if line:
                for ii in range(len(self.xarray_)-1):
                    if colors == None:
                        svgobj.line(self.xarray_[ii], self.yarray_[ii], self.xarray_[ii+1], self.yarray_[ii+1],
                                    'black', gid=f"base-{ii}/circleID")
                    else:
                        svgobj.line(self.xarray_[ii], self.yarray_[ii], self.xarray_[ii+1], self.yarray_[ii+1],
                                    colors[ii], gid=f"base-{ii}/circleID")
            else:
                if pairs:
                    for pair in pairs:
                        if sequence and sequence[pair['from']] == '@' and sequence[pair['to']] == '@': continue
                        svgobj.line(offset_x + self.xarray_[pair['from']], offset_y + self.yarray_[pair['from']],
                         offset_x + self.xarray_[pair['to']], offset_y + self.yarray_[pair['to']],
                          pair['color'], self.NODE_R, alpha=min([alpha[pair['from']],alpha[pair['to']]]), gid=f"pairID-{pair['from']}-{pair['to']}")

                for ii in range(0,len(self.xarray_)):
                    if sequence and sequence[ii] == '@': continue
                    if colors == None:
                        svgobj.circle(self.xarray_[ii] + offset_x, self.yarray_[ii] + offset_y, self.NODE_R, "#000000", "#000000", alpha[ii], gid=f"base-{ii}/centerID")
                    else:
                        svgobj.circle(self.xarray_[ii] + offset_x, self.yarray_[ii] + offset_y, self.NODE_R, colors[ii], colors[ii], alpha[ii], gid=f"base-{ii}/centerID")


                text_offset_x = 0
                text_offset_y = 0
                text_size = self.NODE_R * 1.5
                if svg_mode:
                    text_offset_x = -4.0
                    text_offset_y = (text_size)/2.0 - 1.0
                if sequence and render_in_letter:

                    # write 5' 3' markers
                    if show_primes:
                        text = svgobj.text(self.xarray_[0] + offset_x - math.sin(external_offset)*2.5*self.NODE_R,
                         self.yarray_[0] + offset_y - math.cos(external_offset)* 2.5*self.NODE_R, self.NODE_R * 1.5, "#000000", "center", "5'", 1, gid='5marker')
                        text = svgobj.text(self.xarray_[-1] + offset_x - math.sin(external_offset)*2.5*self.NODE_R,
                         self.yarray_[-1] + offset_y - math.cos(external_offset)* 2.5*self.NODE_R, self.NODE_R * 1.5, "#000000", "center", "3'", 1, gid='3marker')

                    draw_lines = []
                    for ii in range(0,len(self.xarray_)):
                        nucleo_text_size = self.NODE_R * 1.1
                        if sequence[ii] == '@':
                            draw_lines.append(ii)
                            continue
                        if colors[ii] == [0,0,0]:
                            color = "#FFFFFF"
                        else:
                            color = "#000000"
                        if svg_mode:
                            # text_offset_x = -4.0
                            # text_offset_y = (text_size)/2.0 - 1.0
                            svgobj.text(self.xarray_[ii] + offset_x + text_offset_x, self.yarray_[ii] + offset_y + text_offset_y, nucleo_text_size, color, "center", sequence[ii], gid=f"base-{ii}/seqID")
                        else:
                            svgobj.text(self.xarray_[ii] + offset_x, self.yarray_[ii] + offset_y-1, nucleo_text_size, color, "center", sequence[ii], alpha[ii], gid=f"base-{ii}/seqID")
                    # Add sequence numbering
                    text_numbering = []
                    if sequence and (numbering is not None):
                        if len(numbering[0]) != len(sequence):
                            raise RuntimeError("Need to have the same number of nucleotide numbers as sequence letters.")
                        numbering, numbering_spacing = numbering[0], numbering[1]
                        for ii in range(len(numbering)):
                            if sequence[ii] == '@': continue
                            if numbering[ii] % numbering_spacing == 0:
                                [x, y] = self.get_distant_xy_pos(ii, offset_x + text_offset_x, offset_y + text_offset_y, dist_scaler=4.7)
                                # [x, y] = self.get_distant_xy_pos(ii, offset_x , offset_y )
                                text = svgobj.text(x, y, text_size*1.25, "#efa053", "center", str(numbering[ii]), gid=f"numbID-{ii}")
                                # text.set_bbox(dict(boxstyle="round,pad=0.05", fc="white", ec="white"))
                                # text = svgobj.text(x, y, text_size, "#000000", "center", str(numbering[ii]), gid=f"numbID-{ii}")
                                text_numbering.append((text, ii))
                                # svgobj.text(x, y, text_size*.5, "#b9c4c2", "center", str(numbering[ii]))
                    text_custom = []
                    if custom_text:
                        for ii, custom in enumerate(custom_text):
                            [x, y] = self.get_distant_xy_pos(custom[0], offset_x + text_offset_x, offset_y + text_offset_y, dist_scaler=7.5)
                            # text = svgobj.text(x, y, text_size*.8, "#000000", "center", custom[1], gid=f"customID-{ii}", check_position=False)
                            # text = svgobj.text(x, y, text_size*1.4, "#000000", "center", custom[1], gid=f"customID-{ii}", check_position=False)
                            text = svgobj.text(x, y, text_size*1.35, "#000000", "center", custom[1], gid=f"customID-{ii}", check_position=False)
                            text_custom.append((text, custom))
                    check_overlap(svgobj.fig, svgobj.ax)
                    if custom_line:
                        for text, custom in text_custom:
                            if len(custom) < 3: continue
                            new_x, new_y = text.get_position()
                            nucl1 = self.xarray_[custom[2][0]] + offset_x + text_offset_x, self.yarray_[custom[2][0]] + offset_y + text_offset_y
                            nucl2 = self.xarray_[custom[2][1]] + offset_x + text_offset_x, self.yarray_[custom[2][1]] + offset_y + text_offset_y
                            mx = (nucl1[0] + nucl2[0])/2
                            my = (nucl1[1] + nucl2[1])/2
                            svgobj.line(
                                new_x, new_y, mx, my, custom_line_col, width=custom_line_width, alpha=custom_line_alpha,
                                linestyle=custom_line_style
                            )
                        for text, ii in text_numbering:
                            new_x, new_y = text.get_position()
                            nucl = self.xarray_[ii] + offset_x + text_offset_x, self.yarray_[ii] + offset_y + text_offset_y
                            svgobj.line(
                                new_x, new_y, nucl[0], nucl[1], custom_line_col, width=custom_line_width, alpha=custom_line_alpha,
                                linestyle=custom_line_style
                            )
                        added_lines = []
                        for ii in draw_lines:
                            if ii-1 in added_lines: continue
                            added_lines.append(ii)
                            x1, y1 = self.xarray_[ii-1] + offset_x + text_offset_x, self.yarray_[ii-1] + offset_y + text_offset_y
                            x2, y2 = self.xarray_[ii+1] + offset_x + text_offset_x, self.yarray_[ii+1] + offset_y + text_offset_y
                            svgobj.line(
                                x1, y1, x2, y2, '000000', width=custom_line_width, alpha=custom_line_alpha,
                                linestyle=custom_line_style
                            )

                            # svgobj.text(x, y, text_size*.8, "#000000", "center", custom[1], gid=f"customID-{ii}", check_position=False)
                            # if not svgobj.text(x, y, text_size*.8, "#000000", "center", custom[1], gid=f"customID-{ii}", check_position=True):
                                # [x, y] = self.get_distant_xy_pos(custom[0], offset_x + text_offset_x, offset_y + text_offset_y, dist_scaler=5)
                                # svgobj.text(x, y, text_size*.8, "#000000", "center", custom[1], gid=f"customID-{ii}", check_position=True)
                                # pass


    def get_distant_xy_pos(self, idx, offset_x, offset_y, dist_scaler=3.5):
        # neighbor_idx = idx - 1
        # if idx == 0:
        #     neighbor_idx = idx + 1
        # perp_vec = (-self.yarray_[idx] + self.yarray_[neighbor_idx], self.xarray_[idx] - self.xarray_[neighbor_idx])
        # perp_angle = math.atan2(perp_vec[0], perp_vec[1])

        # trial_angles = [perp_angle, -perp_angle]
        trial_angles = [0, np.pi/2, np.pi, 3 * np.pi/2]
        max_min_dist = -1
        max_min_angle = 0

        for angle in trial_angles:
            x_pos = self.xarray_[idx] - math.sin(angle)*dist_scaler*self.NODE_R
            y_pos = self.yarray_[idx] - math.cos(angle)*dist_scaler*self.NODE_R

            min_dist = - 1
            for ii in range(len(self.xarray_)):
                if abs(ii - idx) < 2:
                    continue
                dist = math.sqrt((x_pos - self.xarray_[ii])**2 + (y_pos - self.yarray_[ii])**2)
                if min_dist == -1 or dist < min_dist:
                    min_dist = dist

            if max_min_dist == -1 or min_dist > max_min_dist:
                max_min_dist = min_dist
                max_min_angle = angle

        x_pos = self.xarray_[idx] + offset_x - math.sin(max_min_angle)*dist_scaler*self.NODE_R
        y_pos = self.yarray_[idx] + offset_y - math.cos(max_min_angle)*dist_scaler*self.NODE_R
        return [x_pos, y_pos]

    def get_coords(self, xarray, yarray, PRIMARY_SPACE, PAIR_SPACE):

        if(self.root_ != None) :
            get_coords_recursive(self.root_, xarray, yarray, PRIMARY_SPACE, PAIR_SPACE)
        else :
            for ii in range(0,len(xarray)):
                xarray[ii] = 0
                yarray[ii] = ii * PRIMARY_SPACE

    def setup_coords(self, NODE_R, PRIMARY_SPACE, PAIR_SPACE, external_multiplier, external_offset):
        if self.root_ != None:
            setup_coords_recursive(self.root_, None, 0, 0, 0, 1, NODE_R, PRIMARY_SPACE, PAIR_SPACE, external_multiplier, external_offset)




def check_overlap(fig, ax, max_x_shift = 0.5, max_y_shift = 0.3):
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
        # if text.get_text().upper() in ['A', 'G', 'C', 'U', 'T', "5'", "3'"]:
        if text.get_text().upper() in ['A', 'G', 'C', 'U', 'T']:
            is_moveable = False
        if text.get_text().isdigit():
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
