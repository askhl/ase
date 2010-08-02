# encoding: utf-8

"Module for homogeneous deformation and calculations of elastic constants."

import gtk
from ase.gui.simulation import Simulation
from ase.gui.minimize import MinimizeMixin
from ase.gui.energyforces import OutputFieldMixin
from ase.gui.widgets import oops, pack, AseGuiCancelException
import ase
import numpy as np

class HomogeneousDeformation(Simulation, MinimizeMixin, OutputFieldMixin):
    "Window for homogeneous deformation and elastic constants."
    
    def __init__(self, gui):
        Simulation.__init__(self, gui)
        self.set_title("Homogeneous scaling")
        vbox = gtk.VBox()
        self.packtext(vbox, "XXXX Bla bla bla.")
        self.packimageselection(vbox, txt1="", txt2="")
        self.start_radio_nth.set_active(True)
        pack(vbox, gtk.Label(""))

        # Radio buttons for choosing deformation mode.
        tbl = gtk.Table(4,3)
        for i, l in enumerate(('3D', '2D', '1D')):
            l = l + " deformation   "
            lbl = gtk.Label(l)
            tbl.attach(lbl, i, i+1, 0, 1)
        self.radio_bulk = gtk.RadioButton(None, "Bulk")
        tbl.attach(self.radio_bulk, 0, 1, 1, 2)
        self.radio_xy = gtk.RadioButton(self.radio_bulk, "xy-plane")
        tbl.attach(self.radio_xy, 1, 2, 1, 2)
        self.radio_xz = gtk.RadioButton(self.radio_bulk, "xz-plane")
        tbl.attach(self.radio_xz, 1, 2, 2, 3)
        self.radio_yz = gtk.RadioButton(self.radio_bulk, "yz-plane")
        tbl.attach(self.radio_yz, 1, 2, 3, 4)
        self.radio_x = gtk.RadioButton(self.radio_bulk, "x-axis")
        tbl.attach(self.radio_x, 2, 3, 1, 2)
        self.radio_y = gtk.RadioButton(self.radio_bulk, "y-axis")
        tbl.attach(self.radio_y, 2, 3, 2, 3)
        self.radio_z = gtk.RadioButton(self.radio_bulk, "z-axis")
        tbl.attach(self.radio_z, 2, 3, 3, 4)
        tbl.show_all()
        pack(vbox, [tbl])
        self.deformtable = [
            (self.radio_bulk, (1,1,1)),
            (self.radio_xy, (1,1,0)),
            (self.radio_xz, (1,0,1)),
            (self.radio_yz, (0,1,1)),
            (self.radio_x, (1,0,0)),
            (self.radio_y, (0,1,0)),
            (self.radio_z, (0,0,1))]
        self.deform_label = gtk.Label("")
        pack(vbox, [self.deform_label])
        self.choose_possible_deformations(first=True)

        # Parameters for the deformation
        framedef = gtk.Frame("Deformation:")
        vbox2 = gtk.VBox()
        vbox2.show()
        framedef.add(vbox2)
        self.max_scale = gtk.Adjustment(0.010, 0.001, 10.0, 0.001)
        max_scale_spin = gtk.SpinButton(self.max_scale, 10.0, 3)
        pack(vbox2, [gtk.Label("Maximal scale factor: "), max_scale_spin])
        self.scale_offset = gtk.Adjustment(0.0, -10.0, 10.0, 0.001)
        scale_offset_spin = gtk.SpinButton(self.scale_offset, 10.0, 3)
        pack(vbox2, [gtk.Label("Scale offset: "), scale_offset_spin])
        self.nsteps = gtk.Adjustment(5, 3, 100, 1)
        nsteps_spin = gtk.SpinButton(self.nsteps, 1, 0)
        pack(vbox2, [gtk.Label("Number of steps: "), nsteps_spin])
        
        # Atomic relaxations
        framerel = gtk.Frame("Atomic relaxations:")
        vbox2 = gtk.VBox()
        vbox2.show()
        framerel.add(vbox2)
        self.radio_relax_on = gtk.RadioButton(None, "On   ")
        self.radio_relax_off = gtk.RadioButton(self.radio_relax_on, "Off")
        self.radio_relax_off.set_active(True)
        pack(vbox2, [self.radio_relax_on, self.radio_relax_off])
        self.make_minimize_gui(vbox2)
        for r in (self.radio_relax_on, self.radio_relax_off):
            r.connect("toggled", self.relax_toggled)
        self.relax_toggled()
        pack(vbox, [framedef, gtk.Label(" "), framerel])
        pack(vbox, gtk.Label(""))
        
        # Results
        pack(vbox, [gtk.Label("Results:")])
        self.radio_results_optimal = gtk.RadioButton(
            None, "Load optimal configuration (XXX broken!)")
        self.radio_results_all =  gtk.RadioButton(
            self.radio_results_optimal, "Load all configurations")
        self.radio_results_optimal.set_active(True)
        pack(vbox, [self.radio_results_optimal])
        pack(vbox, [self.radio_results_all])

        # Output field
        outframe = self.makeoutputfield(None)
        fitframe = gtk.Frame("Fit:")
        vbox2 = gtk.VBox()
        vbox2.show()
        fitframe.add(vbox2)
        self.radio_fit_2 = gtk.RadioButton(None, "2nd")
        self.radio_fit_3 = gtk.RadioButton(self.radio_fit_2, "3rd")
        self.radio_fit_2.connect("toggled", self.change_fit)
        self.radio_fit_3.connect("toggled", self.change_fit)
        self.radio_fit_3.set_active(True)
        pack(vbox2, [gtk.Label("Order of fit: "), self.radio_fit_2,
                     self.radio_fit_3])
        pack(vbox2, [gtk.Label("")])
        scrwin = gtk.ScrolledWindow()
        scrwin.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.fit_output = gtk.TextBuffer()
        txtview = gtk.TextView(self.fit_output)
        txtview.set_editable(False)
        scrwin.add(txtview)
        scrwin.show_all()
        self.fit_win = scrwin
        print "SIZE RQ:", scrwin.size_request()
        vbox2.pack_start(scrwin, True, True, 0)
        hbox = gtk.HBox(homogeneous=True)
        for w in [outframe, fitframe]:
            hbox.pack_start(w)
            w.show()
        pack(vbox, hbox)    
        pack(vbox, gtk.Label(""))

        # Status field
        self.status_label = gtk.Label("")
        pack(vbox, [self.status_label])
        
        # Run buttons etc.
        self.makebutbox(vbox)
        vbox.show()
        self.add(vbox)
        self.show()
        self.gui.register_vulnerable(self)

    def choose_possible_deformations(self, first=False):
        """Turn on sensible radio buttons.

        Only radio buttons corresponding to deformations in directions
        with periodic boundary conditions should be turned on.
        """
        if self.setup_atoms():
            pbc = self.atoms.get_pbc()
        else:
            pbc = [False, False, False]
        for radio, requirement in self.deformtable:
            ok = True
            for i in range(3):
                if requirement[i] and not pbc[i]:
                    ok = False
            radio.set_sensitive(ok)
            if first and ok:
                # The first acceptable choice, choose it to prevent
                # inconsistent state.
                radio.set_active(True)
                first = False

    def relax_toggled(self, *args):
        "Turn minimization widgets on or off."
        state = self.radio_relax_on.get_active()
        for widget in (self.algo, self.fmax_spin, self.steps_spin):
            widget.set_sensitive(state)
        
    def notify_atoms_changed(self):
        "When atoms have changed, check for the number of images."
        self.setupimageselection()
        self.choose_possible_deformations()

    def get_deformation_axes(self):
        """Return which axes the user wants to deform along."""
        for but, deform in self.deformtable:
            if but.get_active():
                return np.array(deform)
        # No deformation chosen!
        oops("No deformation chosen: Please choose a deformation mode.")
        return False
            
    def run(self, *args):
        """Make the deformation."""
        self.output.set_text("")
        if not self.setup_atoms():
            return
        deform_axes = self.get_deformation_axes()
        if deform_axes is False:
            return  #Nothing to do!

        # Prepare progress bar
        if self.radio_relax_on.get_active():
            fmax = self.fmax.value
            mininame = self.minimizers[self.algo.get_active()]
            self.begin(mode="scale/min", algo=mininame, fmax=fmax,
                       steps=self.steps.value, scalesteps=self.nsteps.value)
        else:
            self.begin(mode="scale", scalesteps=self.nsteps.value)
        try:
            logger_func = self.gui.simulation['progress'].get_logger_stream
        except (KeyError, AttributeError):
            logger = None
        else:
            logger = logger_func()  # Don't catch errors in the function.

        # Display status message
        self.status_label.set_text("Running ...")
        self.status_label.modify_fg(gtk.STATE_NORMAL,
                                    gtk.gdk.color_parse('#AA0000'))
        while gtk.events_pending():
            gtk.main_iteration()

        # Do the scaling
        scale = self.max_scale.value
        steps = np.linspace(-scale, scale, self.nsteps.value)
        steps += self.scale_offset.value
        undef_cell = self.atoms.get_cell()
        results = []
        txt = "Strain\t\tEnergy [eV]\n"
        # If we load all configurations, prepare it.
        if self.radio_results_all.get_active():
            self.prepare_store_atoms()

        try:
            # Now, do the deformation
            for i, d in enumerate(steps):
                deformation = np.diag(1.0 + d * deform_axes)
                self.atoms.set_cell(np.dot(undef_cell, deformation),
                                    scale_atoms=True)
                if self.gui.simulation.has_key('progress'):
                    self.gui.simulation['progress'].set_scale_progress(i)
                if self.radio_relax_on.get_active():
                    algo = getattr(ase.optimize, mininame)
                    minimizer = algo(self.atoms, logfile=logger)
                    minimizer.run(fmax=fmax, steps=self.steps.value)
                e = self.atoms.get_potential_energy()
                results.append((d, e))
                txt = txt + ("%.3f\t\t%.3f\n" % (d, e))
                self.output.set_text(txt)
                if self.radio_results_all.get_active():
                    self.store_atoms()
        except AseGuiCancelException:
            # Update display to reflect cancellation of simulation.
            self.status_label.set_text("Calculation CANCELLED.")
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#AA4000'))
        except MemoryError:
            self.status_label.set_text("Out of memory, consider using LBFGS instead")
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#AA4000'))
            
        else:
            # Update display to reflect succesful end of simulation.
            self.status_label.set_text("Calculation completed.")
            self.status_label.modify_fg(gtk.STATE_NORMAL,
                                        gtk.gdk.color_parse('#007700'))
                     
        if results:
            self.do_fit(np.array(results))
            if self.radio_results_optimal.get_active():
                if self.minimum_ok:
                    deformation = np.diag(1.0 + self.x0 * deform_axes)
                    self.atoms.set_cell(np.dot(undef_cell, deformation),
                                        scale_atoms=True)
                    if self.radio_relax_on.get_active():
                        if self.gui.simulation.has_key('progress'):
                            self.gui.simulation['progress'].set_scale_progress(
                                len(steps))
                        algo = getattr(ase.optimize, mininame)
                        minimizer = algo(self.atoms, logfile=logger)
                        minimizer.run(fmax=fmax, steps=self.steps.value)
                    # Store the optimal configuration.
                    self.prepare_store_atoms()
                    self.store_atoms()  
                else:
                    oops("No trustworthy minimum: Old configuration kept.")
            self.activate_output()
            self.gui.notify_vulnerable()
        self.end()    
            
        # If we store all configurations: Open movie window and energy graph
        if self.gui.images.nimages > 1:
            self.gui.movie()
            assert not np.isnan(self.gui.images.E[0])
            if not self.gui.plot_graphs_newatoms():
                expr = 'i, e - E[-1]'            
                self.gui.plot_graphs(expr=expr)
            # Continuations should use the best image
            nbest = np.argmin(np.array(results)[:,1])
            self.start_nth_adj.value = nbest
            

    def change_fit(self, widget):
        "Repeat the fitting if the order is changed."
        # It may be called both for the button being turned on and the
        # one being turned off.  But we only want to call do_fit once.
        # And only if there are already cached results (ie. if the
        # order is changed AFTER the calculation is done).
        if widget.get_active() and getattr(self, "results", None) is not None:
            self.do_fit(None)
                
    def do_fit(self, results):
        "Fit the results to a polynomial"
        if results is None:
            results = self.results  # Use cached results
        else:
            self.results = results  # Keep for next time
        self.minimum_ok = False
        if self.radio_fit_3.get_active():
            order = 3
        else:
            order = 2
            
        if len(results) < 3:
            txt = ("Insufficent data for a fit\n(only %i data points)\n"
                   % (len(results),) )
            order = 0
        elif len(results) == 3 and order == 3:
            txt = "REVERTING TO 2ND ORDER FIT\n(only 3 data points)\n\n"
            order = 2
        else:
            txt = ""

        if order > 0:
            fit0 = np.poly1d(np.polyfit(results[:,0], results[:,1], order))
            fit1 = np.polyder(fit0, 1)
            fit2 = np.polyder(fit1, 1)

            x0 = None
            for t in np.roots(fit1):
                if fit2(t) > 0:
                    x0 = t
                    break
            if x0 is None:
                txt = txt + "No minimum found!"
            else:
                e0 = fit0(x0)
                e2 = fit2(x0)
                txt += "E = "
                if order == 3:
                    txt += "A(x - x0)³ + "
                txt += "B(x - x0)² + C\n\n"
                txt += "B = %.5g eV\n" % (e2,)
                txt += "C = %.5g eV\n" % (e0,)
                txt += "x0 = %.5g\n" % (x0,)
                lowest = self.scale_offset.value - self.max_scale.value
                highest = self.scale_offset.value + self.max_scale.value
                if x0 < lowest or x0 > highest:
                    txt += "\nWARNING: Minimum is outside interval\n"
                    txt += "It is UNRELIABLE!\n"
                else:
                    self.minimum_ok = True
                    self.x0 = x0
        self.fit_output.set_text(txt)
        
