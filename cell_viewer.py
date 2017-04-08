#!/usr/bin/env python

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

import inspect
import Tkinter
import PIL
import PIL.ImageTk 
import cairo
import pyradbar
import numpy.random
import random

class CellViewer(object):
    def __init__(self, parent,  cell=None):
        # Setup member variables.
        self._parent = parent
        self._widget = {}
        self._cell = cell
        
        # Setup frames.
        self._cell_frame = Tkinter.Frame(self._parent)
        self._conf_frame = Tkinter.Frame(self._parent)
        
        # Set drawing area height and width.
        self._width = 420
        self._height = 420
        
        # Create cairo context.
        self._cairo_surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, self._width, self._height)
        self._cairo_context = cairo.Context(self._cairo_surface)
        
        # Put context to frame.
        self._image_ref = PIL.ImageTk.PhotoImage(PIL.Image.frombuffer("RGBA", (self._width, self._height), self._cairo_surface.get_data(), "raw", "BGRA", 0, 1))
        self._cell_label = Tkinter.Label(self._cell_frame, image=self._image_ref)
        self._cell_label.pack()
        
        # Setup configuration frame.
        self._widget['cell_label'] = Tkinter.Label( self._conf_frame, anchor='w', text='Cell', width=10)
        self._widget['cell_label'].grid(row=0, column=0)
        self._widget['rho_label'] = Tkinter.Label( self._conf_frame, anchor='w', text='Rho:', width=10)
        self._widget['rho_label'].grid(row=1, column=0)
        self._widget['rho'] = Tkinter.Label( self._conf_frame, anchor='w', text='', width=40)
        self._widget['rho'].grid(row=1, column=1)
        self._widget['lambda_label'] = Tkinter.Label( self._conf_frame, anchor='w', text='Lambda:', width=10)
        self._widget['lambda_label'].grid(row=2, column=0)
        self._widget['lambda'] = Tkinter.Label( self._conf_frame, anchor='w', text='', width=40)
        self._widget['lambda'].grid(row=2, column=1)
        self._widget['genus_label'] = Tkinter.Label( self._conf_frame, anchor='w', text='Genus:', width=10)
        self._widget['genus_label'].grid(row=3, column=0)
        self._widget['genus'] = Tkinter.Label( self._conf_frame, anchor='w', text='', width=40)
        self._widget['genus'].grid(row=3, column=1)
        self._widget['num_punc_label'] = Tkinter.Label( self._conf_frame, anchor='w', text='Num Punct:', width=10)
        self._widget['num_punc_label'].grid(row=4, column=0)
        self._widget['num_punc'] = Tkinter.Label( self._conf_frame, anchor='w', text='', width=40)
        self._widget['num_punc'].grid(row=4, column=1)
        #self._widget['cell_text'] = Tkinter.Text(self._conf_frame, background='white', borderwidth='2p', exportselection=0, foreground='black', height=1, width=100, wrap=Tkinter.NONE)
        #self._widget['cell_text'].grid(row=3, column=0)

        self._cell_frame.grid(row=0, column=0)
        self._conf_frame.grid(row=1, column=0)
        
        # Draw initial cell.
        self.update(cell)
    
    def update(self,  e=None):
        # Pick a new cell if there is nothing more to do with the old cell.
        if self._cell is None or self._cell.degree() <= 0:
            valid_permutation_found = False
            while not valid_permutation_found:                
                k = random.randrange(3, 14)
                l = numpy.random.permutation(k)
                self._cell = pyradbar.Cell( l, suspend=False )
                if self._cell.num_fixed_pts_lam() == 0:
                    valid_permutation_found = True
            self.update_wrt_cell()
        else:
            # Take a face of the cell.
            i = random.randrange(0, self._cell.degree() + 1)
            self._cell.face(i)
            self.update_wrt_cell()
    
    def update_wrt_cell(self):
        # Draw the cell and update the configuration widgets.
        pyradbar.draw_cell(self._cairo_context, self._cell )
        self._image_ref = PIL.ImageTk.PhotoImage(PIL.Image.frombuffer("RGBA", (self._width, self._height), self._cairo_surface.get_data(), "raw", "BGRA", 0, 1))
        self._cell_label.configure(image=self._image_ref)
        self._widget['rho'].configure(text=str(self._cell.cyc_rho()))
        self._widget['lambda'].configure(text=str(self._cell.cyc_lam()))
        self._widget['genus'].configure(text=str(self._cell.get_genus()))
        self._widget['num_punc'].configure(text=str(self._cell.get_num_punc()))
        #self._widget['cell_text'].delete('0.0', Tkinter.END)
        #self._widget['cell_text'].insert('0.0', str(self._cell))
    
if __name__ == "__main__":
    # Get the path of the script.
    path = None
    try:
        filename = inspect.getframeinfo(inspect.currentframe()).filename
        path = os.path.dirname(os.path.abspath(filename))
    except:
        path = '.'

    # Setup the Tk window.
    root = Tkinter.Tk()
    root.resizable(width=False, height=False)
    root.wm_title('Cell Viewer')
    try:
        img = Tkinter.PhotoImage(file=path + '/data/icon.png')
        root.tk.call('wm', 'iconphoto', root._w, img)
    except:
        pass
    
    # Create and setup the context of the Tk window.
    cw = CellViewer(root)
    root.bind("<Return>", cw.update)
    
    # Start the mainloop.
    root.mainloop()
