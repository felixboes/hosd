# The software pyradbar is a bunch of programs to compute the homology of
# Sullivan diagrams.
# Copyright (C) 2015 - 2017  Felix Boes
#
# This file is part of pyradbar.
#
# pyradbar is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyradbar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyradbar.  If not, see <http://www.gnu.org/licenses/>.

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk
import cairo
import math
import inspect
import sys
import perm
from cell_colors import cell_color
from cell import *

def centroid_polyhedron(vertices):
    """
    Computes the centroid of a polyhedron.
    Is is necessary to provide the vertices in the correct order, i.e.
    the series of coordinates (x_i, y_i) describes the convex hull.
    
    Compare http://geomalgorithms.com/a01-_area.html#2D%20Polygons
    """
    
    n = len(vertices)
    
    # exclude degenerate cases
    if n < 3:
        if n == 2:
            return 0.5*(vertices[0][0] + vertices[1][0]),  0.5*(vertices[0][1] + vertices[1][1])
        elif n == 1:
            return vertices[0][0],  vertices[0][1]
        else:
            return 0.0,  0.0
    
    # c = [ x_n*y_0 - x_0y_n, x_0*y_1 - x_1*y_0, ... ]
    c = []
    for i in xrange(-1, n - 1):
        c += [ vertices[i][0]*vertices[i+1][1] - vertices[i+1][0]*vertices[i][1] ]
    
    # A = 6 * signed area
    A = 0.0
    for a in c:
        A += a
    A = 3*A
    
    # x,y = coordinates of the centroid
    x = 0.0
    y = 0.0
    for i in xrange(-1, n - 1):
        x += (vertices[i][0] + vertices[i+1][0]) * c[i+1]
        y += (vertices[i][1] + vertices[i+1][1]) * c[i+1]
    
    return x/A,  y/A

def centroid_of_polyhedron_inside_circle( angles ):
    """
    Computes the centroid of the polyhedron spanned by the
    points on the standard circle given by a list of angles.
    
    Returns the coordinates x,y of the centroid.
    """
    coordinates = [ (math.cos(2*math.pi*phi),  math.sin(2*math.pi*phi)) for phi in sorted(angles) ]
    return centroid_polyhedron(coordinates)

def normalized_centroid_cycle(cyclic_perm, num_symb):
    """
    Computes the normaized centroid of a cyclic permutation.
    """
    n = float(num_symb)
    angles = [ x/n for x in cyclic_perm ]
    
    x, y = centroid_of_polyhedron_inside_circle( angles )
    radius = math.sqrt( x**2 + y**2 )
    bump_factor = 0.5
    
    return bump_factor*x,  bump_factor*y

def draw_cell(cr, cell):
    p=cell.degree()+1
    cycles = cell.cyc_lam()
    surface = cell.get_surfaces()
    genus = {}
    punct = {}
    surf_of_cyc = {}
    for i, surf in enumerate(surface):
        for cyc_idx in surf[0]:
            surf_of_cyc[cyc_idx] = i
            genus[cyc_idx] = surf[1]
            punct[cyc_idx] = surf[2]
    
    cr.set_source_rgb(1, 1, 1)
    cr.rectangle(0, 0, 420, 420)
    cr.fill()
    
    cr.set_source_rgb(0, 0, 0)
    cr.set_line_width(2)
    
    radius = 200
    mid = (10 + radius,  10 + radius)
    
    i = 0
    for x in cycles:
        if len(x) == 1:
            draw_degenerate_loop( cr, mid, radius, x[0] / float(p), genus[i], punct[i], surf_of_cyc[i] )
        elif len(x) == 2:
            draw_line( cr, mid, radius, x[0] / float( p ), x[1] / float( p ), genus[i], punct[i], surf_of_cyc[i] )
        elif len(x) > 2:
            draw_cycle( cr, mid, radius, x, p, genus[i], punct[i], surf_of_cyc[i] )
        i += 1
    
    #draw circle
    cr.set_source_rgb(0, 0, 0)
    cr.set_line_width(2)
    cr.arc(mid[0], mid[1], radius, 0,  2*math.pi)
    cr.stroke()
    
def draw_line(cr, mid, radius, phi, psi, genus, punct, line_color_id):
    # setup colors
    white = cell_color('white')
    line_color = cell_color(line_color_id)
    thick_line = [10]
    slim_line  = [ 2]
    back_color = None
    draw_color = None
    
    # set bezier points
    p_0 = mid[0] + radius*math.cos(2*math.pi*phi), mid[1] + radius*math.sin(2*math.pi*phi)
    p_3 = mid[0] + radius*math.cos(2*math.pi*psi), mid[1] + radius*math.sin(2*math.pi*psi)
    t = 2.5 * min( math.fabs( phi - psi),  1 - math.fabs( phi - psi) )
    p_1 = tuple( p_0[i] + t*( mid[i] - p_0[i]) for i in range( 2 ) )
    p_2 = tuple( p_3[i] + t*( mid[i] - p_3[i]) for i in range( 2 ) )
    
    # draw lines
    back_color = white + thick_line
    draw_color = line_color + slim_line
    
    for val in (back_color, draw_color):
        cr.set_source_rgb(val[0], val[1], val[2])
        cr.set_line_width(val[3])
        cr.move_to( *p_0 )
        cr.curve_to( *( list(p_1) + list(p_2) + list(p_3) ) )
        
        cr.stroke()
    
    # compute midpoint of the line and draw weight
    # compare http://stackoverflow.com/questions/6244272/drawing-sections-of-a-bezier-curve
    x, y = p_0[0]/8.0 + 3.0*p_1[0]/8.0 + p_3[0]/8.0 + 3.0*p_2[0]/8.0 , p_0[1]/8.0 + 3.0*p_1[1]/8.0 + p_3[1]/8.0 + 3.0*p_2[1]/8.0
    draw_weight(cr, x, y, genus, punct, radius, line_color_id)

def draw_cycle(cr, mid, radius, cyclic_perm, num_symb, genus, punct, line_color_id):
    # setup colors
    white = cell_color('white')
    line_color = cell_color(line_color_id)
    thick_line = [10]
    slim_line  = [ 2]
    back_color = None
    draw_color = None
    
    #center of the cyclic permutation.
    x, y = normalized_centroid_cycle(cyclic_perm, num_symb)
    
    # compute correct position i.e. scaling and offset
    x = radius * x + mid[0]
    y = radius * y + mid[1]
    
    inner_angle = min( cyclic_perm ) / float(num_symb)
    
    # draw lines
    back_color = white + thick_line
    draw_color = line_color + slim_line
    for phi in cyclic_perm :
        phi = float(phi) / num_symb
        for val in (back_color, draw_color):
            cr.set_source_rgb(val[0], val[1], val[2])
            cr.set_line_width(val[3])
            p_0 = x, y
            p_3 = mid[0] + radius*math.cos(2*math.pi*phi), mid[1] + radius*math.sin(2*math.pi*phi)
            # t is a bump factor depending on the distance of the normalized centroid of the cycle to the mid point.
            #t = ( 1- math.sqrt( (mid[0] - x)**2 + (mid[1] - y)**2 )/radius )
            t=max( 0.5, len(cyclic_perm) / float(num_symb) )
            p_1 = x + 0.8*t*radius*math.cos(2*math.pi*inner_angle), y + 0.8*t*radius*math.sin(2*math.pi*inner_angle)
            p_2 = tuple( p_3[i] + t*( mid[i] - p_3[i]) for i in range( 2 ) )
            
            cr.move_to( *p_0 )
            cr.curve_to( *( list(p_1) + list(p_2) + list(p_3) ) )
            
            cr.stroke()
        inner_angle = inner_angle + 1.0/len(cyclic_perm)
    
    # draw circle
    draw_weight(cr, x, y, genus, punct, radius, line_color_id)
    
def draw_degenerate_loop(cr, mid, radius, phi, genus, punct, line_color_id):
    # setup colors
    white = cell_color('white')
    line_color = cell_color(line_color_id)
    thick_line = [10]
    slim_line  = [ 2]
    back_color = None
    draw_color = None
    
    #draw line
    back_color = white + thick_line
    draw_color = line_color + slim_line
    for val in (back_color, draw_color):
        cr.set_source_rgb(val[0], val[1], val[2])
        cr.set_line_width(val[3])
        p_0 = mid[0] +     radius*math.cos(2*math.pi*phi), mid[1] +     radius*math.sin(2*math.pi*phi)
        p_1 = mid[0] + 0.9*radius*math.cos(2*math.pi*phi), mid[1] + 0.9*radius*math.sin(2*math.pi*phi)
        
        cr.move_to( *p_0 )
        cr.line_to( *p_1 )
        cr.stroke()
    
    # draw weight
    x, y = mid[0] + 0.85*radius*math.cos(2*math.pi*phi), mid[1] + 0.85*radius*math.sin(2*math.pi*phi)
    draw_weight(cr, x, y, genus, punct, radius, line_color_id)

def draw_weight(cr, x, y, genus, punct, outer_radius, line_color_id):
    white = cell_color('white')
    line_color = cell_color(line_color_id)
    thick_line = [10]
    slim_line  = [ 2]
    back_color = None
    draw_color = None
    
    #draw circle
    back_color = line_color + thick_line
    draw_color = white + slim_line
    for val in (back_color, draw_color):
        cr.set_source_rgb(val[0], val[1], val[2])
        r = (0.05 if val[3] == thick_line[0] else 0.04)*outer_radius
        
        cr.arc(x, y, r, 0, 2*math.pi)
        cr.fill()
        cr.stroke()
    
    #write weight
    draw_color = line_color
    
    cr.set_source_rgb(*draw_color)
    cr.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    cr.set_font_size(8)
    #center text at x,y
    text = str(genus) + ';' + str(punct)
    x_bearing, y_bearing, ext_width, ext_height, x_adv, y_adv = cr.text_extents(text)
    x -= (ext_width/2  + x_bearing)
    y -= (ext_height/2 + y_bearing)
    cr.move_to(x, y)
    cr.show_text(text)
    cr.stroke()

class DrawCell( object ):
    def __init__(self, cell):
        self.cell = cell
        
        surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, 420, 420)
        cr = cairo.Context(surface)
        
        draw_cell(cr, self.cell)
        try:
            surface.write_to_png ('./results/g_' + str(self.cell.get_genus()) + '_m_' + str(self.cell.get_num_punc()) + '/' + self.cell.get_filename_str() +'.png')
        except:
            # Print the error.
            frameinfo = inspect.getframeinfo(inspect.currentframe())
            print frameinfo.filename, frameinfo.lineno
            e, p, t = sys.exc_info()
            print e, p
