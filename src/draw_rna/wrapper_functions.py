from draw_rna.ipynb_draw import draw_struct
import os
import re
from RNA import RNA
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

def get_default_settings():
    return {
        'predict':False,
        'large_mode': False,
        'black_white': False,
        'store_plt_svg': True,
        'file_name': None,
        'numbering': False,
        'segment_lenghts': True,
        'show_energies': False,
        'energies_individual': False,
        'custom_line': False,
        'constraints': None,
        'line': False,
        'add_fake': [],
    }
def create_figure(file_path, seq=None, db=None, colors=None,  black_white=False, file_name=None, numbering=False, segment_lenghts=True, custom_text=None,
                  show_energies=False, energies_individual=False, constraints=None, add_fake=None, title=None, energy_title=True, overwrite_colors=None,
                  color_map=None, return_ax=False, **kwargs):
    if custom_text is None: custom_text = []
    if file_path is not None: final_list, final_str_list, seq, db, structures = ct_to_paired_list_count_all(file_path, constraints, **kwargs)
    else: final_list, final_str_list, seq, db, structures = ct_to_paired_list_count_all(None, seq, db, constraints, **kwargs)
    try:
        # free_bases_str = f"free bases, {structures[-1][0]} nm: {structures[-1][0] * 0.59}"
        free_bases_str = f"{structures[-1][0] * 0.59:.2f}"
    except TypeError:
        free_bases_str = f"free bases, {structures[-1][0]}"
    energies, total_energy = get_energy_list(seq, db)
    if energy_title and title is not None:
        title = (title[0]+f"\nenergy: {total_energy}", title[1], title[2])
    if isinstance(colors, str):
        colors = [colors for _ in seq]
    # print(energies)
    # print(total_energy, 'total energy')
    scalar_map = None
    if show_energies and len(energies):
        colors, scalar_map, segment_energy = map_energy_to_color(energies, structures, seq, energies_individual, set_v=None)
        new_energies = []
        for idx, (segment_str, seg_energy) in enumerate(zip(final_str_list, segment_energy)):
            seg2 = segment_str[2]
            if add_fake is not None:
                for fake in add_fake:
                    if fake[0] <= seg2[0]:
                        seg2 = (seg2[0]+fake[1], seg2[1])
                    if fake[0] <= seg2[1]:
                        seg2 = (seg2[0], (seg2[1]+fake[1]))
            final_str_list[idx] = (segment_str[0], segment_str[1]+f"\n{round(seg_energy,2)}∆", seg2)
        # print(final_str_list)
    if constraints:
        if isinstance(constraints, str):
            seq_arr = list(seq)
            for idx, val in enumerate(constraints):
                if val == 'X':
                    seq_arr[idx] = val
            seq = ''.join(seq_arr)
        else:
            for constraint in constraints:
                seq = seq[:constraint] + 'x' + seq[constraint+1:]
    if file_name is None:
        if file_path is None:
            file_name = 'unknown'
        else:
            file_name=file_path.split('.')[0]
    if not segment_lenghts:
        final_str_list = []
    if numbering: numbering=[list(range(1,len(seq)+1)), numbering]
    else: numbering = None
    if black_white: seq = None
    else: seq = seq.upper()
    if overwrite_colors:
        color = overwrite_colors[-1]
        for idx in overwrite_colors[0]:
            colors[idx] = color

    if add_fake is not None:
        for fake in add_fake[::-1]:
            if seq: seq = seq[0:fake[0]]+'@'*fake[1]+seq[fake[0]:]
            if colors: colors = colors[0:fake[0]]+['@']*fake[1]+colors[fake[0]:]
            if numbering: numbering[0] = numbering[0][0:fake[0]]+['@']*fake[1]+numbering[0][fake[0]:]
            if len(fake) == 2:
                db = db[0:fake[0]]+'@'*fake[1]+db[fake[0]:]
            else:
                db = db[0:fake[0]]+fake[2]*fake[1]+db[fake[0]:]
    ax = draw_struct(seq, db, c=colors, cmap='RdBu_r', scalar_map=scalar_map,  custom_text=final_str_list+custom_text, filename=file_name, color_map=color_map,
                numbering=numbering, title=title, return_ax=return_ax, **kwargs)
    if return_ax:
        return ax
    return free_bases_str, energies, final_str_list

def get_ct(file_path):
    with open(file_path, 'r') as file:
        line = file.readline().strip()
        while line.startswith('>'):
            line = file.readline().strip()
        seq = line
        struct = file.readline().strip()
    return seq, struct

def get_energy_list(seq, struct, store_file=None):
    energy_list = {}
    external_values = 0
    total_energy = 0
    with open('tmp.txt', 'w') as file:
        RNA.eval_structure_simple(seq, struct, 1, file)
    if store_file is not None:
        with open(store_file, 'w') as file:
            RNA.eval_structure_simple(seq, struct, 1, file)
    with open('tmp.txt', 'r') as file:
        for line in file.readlines():
            stripped = line.strip()
            val = int(stripped.split(':')[-1])
            total_energy += val
            digits = re.findall(r'\d+', stripped)
            digits = digits[:-1]
            digits = digits[0:2]
            if len(digits) == 0: external_values += val/100
            else:
                digit = int(digits[0])
                energy_list[digit-1] = val/100
                digit = int(digits[1])
                energy_list[digit-1] = val/100
    os.remove('tmp.txt')
    return energy_list, total_energy/100

def print_energy_list(seq, struct):
    energy_list = {}
    external_values = 0
    total_energy = 0
    with open('tmp.txt', 'w') as file:
        RNA.eval_structure_simple(seq, struct, 1, file)
    with open('tmp.txt', 'r') as file:
        for line in file.readlines():
            # print(line.strip())
            stripped = line.strip()
            val = int(stripped.split(':')[-1])
            total_energy += val
            digits = re.findall(r'\d+', stripped)
            digits = digits[:-1]
            digits = digits[0:2]
            # print(digits)
            if len(digits) == 0: external_values += val/100
            else:
                digit = int(digits[0])
                energy_list[digit-1] = val/100
                digit = int(digits[1])
                energy_list[digit-1] = val/100
    print(round(total_energy/100, 2))
    os.remove('tmp.txt')

def ct_to_paired_list_count_all(file_path, seq=None, db=None, constraints=None, predict=False, **kwargs):
    if file_path is not None:
        if predict:
            seq, db = get_ct(file_path)
            db = predict_structure(seq, constraints, **kwargs)
        else:
            seq, db = get_ct(file_path)
    elif seq is not None and predict:
        db = predict_structure(seq, constraints, **kwargs)

    i, j = 0, len(db)-1
    stack = []
    structures = []
    depth = 0
    depth_total = 0
    while i < len(db) and db[i] == '.':
        stack.append((i, db[i]))
        i+=1
    while i < len(db):
        while i < len(db) and (db[i] == '(' or db[i] == '.'):
            stack.append((i, db[i]))
            i += 1
        tmp_struct = []
        cnt = 0
        while len(stack) and stack[-1][1] == '.':
            stack.pop()
            cnt += 1
        while i < len(db) and len(stack) and \
                ((stack[-1][1] == '(' and db[i] == ')') or
                 (len(stack) >= 2 and stack[-1][1] == '.' and stack[-2][1] == '(' and db[i] == ')') or
                 (i+1 < len(db) and stack[-1][1] == '(' and db[i] == '.' and db[i+1] == ')')):
            depth_total += 1
            match = stack.pop()
            if match[1] == '.':
                match = stack.pop()
                cnt+=1
            if db[i] == '.':
                i+=1
                cnt += 1
            tmp_struct.append((match[0], i))
            i+=1
        tmp_struct.append(cnt)  # number of unpaired in center
        if tmp_struct:
            structures.append(tmp_struct)
            depth -= 1
    while len(stack):
        val = stack.pop()
        if val[1] != '.':
            print("Error", val, i)
        if len(structures) > 0  and isinstance(structures[-1][0], int):
            structures[-1][0] += 1
        else:
            structures.append([1])
    final_list = []
    for struc in structures:
        if len(struc) == 1:
            final_list.append(struc[0]+len(stack)-1)
            continue
        unpaired_middle, start, end = struc[-1], struc[0], struc[-2]
        bases = (len(struc)-1)*2+struc[-1]
        length = round(bases*0.59,2)
        middle_bases = struc[int((len(struc)-1)/2)]
        position = int((start[0]-end[0])/2)+end[0]
        final_list.append((position, bases, length, middle_bases))
    str_list = []
    for item in final_list:
        if not isinstance(item, tuple): continue
        # str_list.append((item[0], f"{item[0]}\n{item[1]}\n{item[2]}", item[3])) # nm bp and middle
        # str_list.append((item[0], f"{item[1]}\n{item[2]}nm", item[3])) # nm and bp
        str_list.append((item[0], f"{item[2]}nm", item[3])) # only nm
    return final_list, str_list, seq, db, structures

def map_energy_to_color(energies, structures, seq, energies_individual=False, set_v=None):
    grey = (169,169,169)
    cmap = plt.get_cmap('seismic')
    energy_vals = np.array(list(energies.values()))
    max_energy = np.max(np.absolute(energy_vals))
    if set_v is None:
        vmin = -max_energy
        vmax = max_energy
    else:
        vmin = set_v[0]
        vmax = set_v[1]


    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    scalarMap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    colors = []
    for char in seq:
        colors.append(grey)
    if energies_individual:
        for segment in structures:
            for pairs in segment[:-1]:
                if pairs[0] in energies:
                    energy = energies[pairs[0]]
                else:
                    energy = energies[pairs[1]]
                colors[pairs[0]] = scalarMap.to_rgba(energy, bytes=True)[:-1]
                colors[pairs[1]] = scalarMap.to_rgba(energy, bytes=True)[:-1]
    else:
        for segment in structures:
            total_energy = 0
            for pairs in segment[:-1]:
                if pairs[0] in energies:
                    energy = energies[pairs[0]]
                else:
                    energy = energies[pairs[1]]
                total_energy += energy
            if total_energy:
                total_energy /= len(segment[:-1])
            for pairs in segment[:-1]:
                colors[pairs[0]] = scalarMap.to_rgba(total_energy, bytes=True)[:-1]
                colors[pairs[1]] = scalarMap.to_rgba(total_energy, bytes=True)[:-1]
    segment_energy = get_energy_per_segment(energies, structures)

    return colors, scalarMap, segment_energy

def get_energy_per_segment(energies, structures):
    segment_energy = []
    for segment in structures:
        total_energy = 0
        for pairs in segment[:-1]:
            if pairs[0] in energies:
                energy = energies[pairs[0]]
            else:
                energy = energies[pairs[1]]
            total_energy += energy
        if total_energy:
            # total_energy /= len(segment[:-1])
            segment_energy.append(total_energy)
    return segment_energy

def predict_structure(sequence, constraints=None, salt=None, **kwargs):
    db_input = '.'*len(sequence)
    if constraints and not isinstance(constraints, str):
        db_input = list(db_input)
        for constraint in constraints:
            db_input[constraint] = 'x'
        db_input = ''.join(db_input)
    elif constraints:
        print(constraints)
        db_input = constraints
    md = RNA.md()
    if salt:
        md.salt = salt
    fc = RNA.fold_compound(sequence, md)
    # fc.hc_add_from_db(db_input) # Stopped working new fix
    for idx, val in enumerate(list(db_input)):
        if val == 'X':
            fc.hc_add_up(idx+1)
    (ss, mfe) = fc.mfe()
    print("energy", mfe)
    return ss



def scale_figure_with_text(fig, scale_factor):
    """Scale a matplotlib figure including all text elements."""
    # Scale figure
    w, h = fig.get_size_inches()
    fig.set_size_inches(w * scale_factor, h * scale_factor)

    # Scale text for all axes
    for ax in fig.get_axes():
        # Title and axis labels
        items = [ax.title, ax.xaxis.label, ax.yaxis.label]

        # Tick labels
        items.extend(ax.get_xticklabels())
        items.extend(ax.get_yticklabels())

        for line in ax.get_lines():
            line.set_linewidth(line.get_linewidth() * scale_factor)
            if line.get_markersize() is not None:
                line.set_markersize(line.get_markersize() * scale_factor)
        # Legend
        if ax.get_legend() is not None:
            items.extend(ax.get_legend().get_texts())
        # Annotations and text elements
        for child in ax.get_children():
            if isinstance(child, plt.Text):
                items.append(child)

        # Scale font sizes
        for item in items:
            item.set_fontsize(item.get_fontsize() * scale_factor)

def get_structure_depth(structures):
    structure_depth = []
    for s in structures[:-1]:
        x1,x2 = s[-2][0], s[-2][1]
        depth = 0
        for s2 in structures[:-1]:
            if s == s2: continue
            y1,y2 = s2[-2][0], s2[-2][1]
            if y1 < x1 and y2 > x2:
                depth += 1
        structure_depth.append(depth)
    return structure_depth