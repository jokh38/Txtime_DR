from raystation import *
from raystation import get_current, await_user_input, CompositeAction, set_progress
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import math, random
import numpy as np

class GUIProgram:
    def __init__(self, master):
        self.orig_exam_name = []
        self.s4DGlist = []
        self.list_mcPN = []
        self.dr_results = {}
        self.gating_params = {'T_on': 0, 'T_off': 0, 'efficiency': 0}
        
        self.master = master
        self.master.title("Total treatment time calculation with dynamic range")
        self.create_widgets()
        self.update_gating_calculation()

    def create_widgets(self):
        case = get_current("Case")
        
        # Frame: Plan Selection
        self.frame_plan_select = ttk.LabelFrame(self.master, text="Plan Selection")
        self.frame_plan_select.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        self.var_plan_name = tk.StringVar(value=case.TreatmentPlans[0].Name)
        options_plan = [i_plan.Name for i_plan in case.TreatmentPlans]

        plan_label = ttk.Label(self.frame_plan_select, text="[PLAN]")
        plan_label.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
        
        plan_dropdown = ttk.OptionMenu(self.frame_plan_select, self.var_plan_name, options_plan[0], *options_plan)
        plan_dropdown.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")

        # Frame: Respiratory Parameters
        self.frame_resp_params = ttk.LabelFrame(self.master, text="Respiratory Parameters")
        self.frame_resp_params.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")

        self.var_resp_period = tk.StringVar(value='4')
        resp_period_label = ttk.Label(self.frame_resp_params, text="[Resp. period (sec)]")
        resp_period_label.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
        resp_period_entry = tk.Entry(self.frame_resp_params, textvariable=self.var_resp_period)
        resp_period_entry.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")
        resp_period_entry.bind('<KeyRelease>', self.update_gating_calculation)

        self.var_layer_switch_time = tk.StringVar(value='2')
        layer_switch_label = ttk.Label(self.frame_resp_params, text="[Layer switch time (sec)]")
        layer_switch_label.grid(row=2, column=0, padx=5, pady=5, sticky="nsew")
        layer_switch_entry = tk.Entry(self.frame_resp_params, textvariable=self.var_layer_switch_time)
        layer_switch_entry.grid(row=3, column=0, padx=5, pady=5, sticky="nsew")

        # Frame: DR Analysis Range
        self.frame_dr_range = ttk.LabelFrame(self.master, text="DR Analysis Range")
        self.frame_dr_range.grid(row=0, column=2, padx=10, pady=10, sticky="nsew")

        self.var_dr_start = tk.StringVar(value='10')
        dr_start_label = ttk.Label(self.frame_dr_range, text="[DR Start]")
        dr_start_label.grid(row=0, column=0, padx=5, pady=2, sticky="nsew")
        dr_start_entry = tk.Entry(self.frame_dr_range, textvariable=self.var_dr_start)
        dr_start_entry.grid(row=1, column=0, padx=5, pady=2, sticky="nsew")

        self.var_dr_end = tk.StringVar(value='200')
        dr_end_label = ttk.Label(self.frame_dr_range, text="[DR End]")
        dr_end_label.grid(row=2, column=0, padx=5, pady=2, sticky="nsew")
        dr_end_entry = tk.Entry(self.frame_dr_range, textvariable=self.var_dr_end)
        dr_end_entry.grid(row=3, column=0, padx=5, pady=2, sticky="nsew")

        self.var_dr_step = tk.StringVar(value='20')
        dr_step_label = ttk.Label(self.frame_dr_range, text="[DR Step]")
        dr_step_label.grid(row=4, column=0, padx=5, pady=2, sticky="nsew")
        dr_step_entry = tk.Entry(self.frame_dr_range, textvariable=self.var_dr_step)
        dr_step_entry.grid(row=5, column=0, padx=5, pady=2, sticky="nsew")

        # Frame: Analysis Control
        self.frame_analysis_control = ttk.LabelFrame(self.master, text="Analysis Control")
        self.frame_analysis_control.grid(row=0, column=3, padx=10, pady=10, sticky="nsew")

        analyze_btn = ttk.Button(self.frame_analysis_control, text="Analyze DR Range", 
                                command=self.analyze_dr_range)
        analyze_btn.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        export_btn = ttk.Button(self.frame_analysis_control, text="Export Results", 
                               command=self.export_results)
        export_btn.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")

        # Frame: Phase vs Amplitude
        self.frame_phase_amp = ttk.LabelFrame(self.master, text="Phase vs Amplitude")
        self.frame_phase_amp.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")
        
        self.var_phase_amp = tk.IntVar(value=1)
        ttk.Radiobutton(self.frame_phase_amp, text="Phase", variable=self.var_phase_amp, 
                       value=0, command=self.update_gating_calculation).grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(self.frame_phase_amp, text="Amplitude", variable=self.var_phase_amp, 
                       value=1, command=self.update_gating_calculation).grid(row=1, column=0, sticky="w")
        
        # Frame: Respiratory Signal
        self.frame_resp_signal = ttk.LabelFrame(self.master, text="Respiratory Signal")
        self.frame_resp_signal.grid(row=1, column=1, padx=10, pady=10, sticky="nsew")
        
        self.var_resp_pattern = tk.IntVar(value=3)
        ttk.Radiobutton(self.frame_resp_signal, text="Sine", variable=self.var_resp_pattern, 
                       value=1, command=self.update_resp_signal).grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(self.frame_resp_signal, text="Tail", variable=self.var_resp_pattern, 
                       value=2, command=self.update_resp_signal).grid(row=1, column=0, sticky="w")
        ttk.Radiobutton(self.frame_resp_signal, text="Long-tail", variable=self.var_resp_pattern, 
                       value=3, command=self.update_resp_signal).grid(row=2, column=0, sticky="w")
        
        # Frame: Gating Window
        self.frame_gating_window = ttk.LabelFrame(self.master, text="Gating Window")
        self.frame_gating_window.grid(row=1, column=2, padx=10, pady=10, sticky="nsew")
        
        self.var_gating_window = tk.IntVar(value=20)
        ttk.Radiobutton(self.frame_gating_window, text="5%", variable=self.var_gating_window, 
                       value=5, command=self.update_gating_calculation).grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(self.frame_gating_window, text="10%", variable=self.var_gating_window, 
                       value=10, command=self.update_gating_calculation).grid(row=1, column=0, sticky="w")
        ttk.Radiobutton(self.frame_gating_window, text="15%", variable=self.var_gating_window, 
                       value=15, command=self.update_gating_calculation).grid(row=2, column=0, sticky="w")
        ttk.Radiobutton(self.frame_gating_window, text="20%", variable=self.var_gating_window, 
                       value=20, command=self.update_gating_calculation).grid(row=3, column=0, sticky="w")

        # Frame: Gating Calculation Results
        self.frame_gating_calc = ttk.LabelFrame(self.master, text="Gating Calculation")
        self.frame_gating_calc.grid(row=1, column=3, padx=10, pady=10, sticky="nsew")

        self.label_ton = ttk.Label(self.frame_gating_calc, text="T_on: 0.00 sec")
        self.label_ton.grid(row=0, column=0, padx=5, pady=2, sticky="w")

        self.label_toff = ttk.Label(self.frame_gating_calc, text="T_off: 0.00 sec")
        self.label_toff.grid(row=1, column=0, padx=5, pady=2, sticky="w")

        self.label_efficiency = ttk.Label(self.frame_gating_calc, text="Efficiency: 0.00")
        self.label_efficiency.grid(row=2, column=0, padx=5, pady=2, sticky="w")

        # Frame: Progress
        self.frame_progress = ttk.LabelFrame(self.master, text="Progress")
        self.frame_progress.grid(row=2, column=0, columnspan=4, padx=10, pady=10, sticky="nsew")

        self.status_text = tk.Label(self.frame_progress, text="Ready for analysis.")
        self.status_text.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        # Frame: Signal Drawing
        self.frame_signal_draw = ttk.LabelFrame(self.master, text="Signal Drawing")
        self.frame_signal_draw.grid(row=3, column=0, columnspan=4, padx=10, pady=10, sticky="nsew")
        
        self.canvas_signal = tk.Canvas(self.frame_signal_draw, width=500, height=150, bg="white")
        self.canvas_signal.pack()

        # Frame: Results Graph
        self.frame_results_graph = ttk.LabelFrame(self.master, text="DR Analysis Results")
        self.frame_results_graph.grid(row=0, column=4, rowspan=5, padx=10, pady=10, sticky="nsew")

        self.canvas_graph = tk.Canvas(self.frame_results_graph, width=800, height=400, bg="white")
        self.canvas_graph.pack(fill=tk.BOTH, expand=True)

        self.update_resp_signal()

    def resp_signal(self, resp_x):
        """Generate respiratory signal patterns based on selected pattern type"""
        if self.var_resp_pattern.get() == 1:
            # Simple sine wave pattern
            resp_y = [math.sin(math.radians(i)) ** 2 for i in resp_x]
        elif self.var_resp_pattern.get() == 2:
            # Tail pattern - asymmetric respiratory curve
            resp_y = [(math.sin(math.radians(i)) + math.sin(math.radians(i) + 3.14) ** 2 + 
                      0.5 * math.sin(math.radians(i) + 1.57) ** 2 + 
                      0.5 * math.sin(math.radians(i) + 0.785) ** 2) for i in resp_x]
        elif self.var_resp_pattern.get() == 3:
            # Long-tail pattern - extended expiration phase
            resp_y = [(math.sin(math.radians(i)) + math.sin(math.radians(i) + 0.223) ** 2 + 
                      0.5 * math.sin(math.radians(i) - 1.13) ** 2) for i in resp_x]
    
        # Normalize to 0-100 range
        arr = np.array(resp_y)
        res = (100 * arr // max(resp_y)).tolist()
        return res

    def update_resp_signal(self):
        """Update respiratory signal visualization on canvas"""
        resp_x = list(range(-400, 550, 1))
        resp_y = self.resp_signal(resp_x)
    
        # Convert to canvas coordinates for display
        tmp_img_sf = 140
        np_rs_y = np.array(resp_y)
        disp_rs_y = np.subtract(tmp_img_sf, np_rs_y).tolist()
    
        # Create line points for canvas drawing
        resp_points = [coord for pair in zip(resp_x, disp_rs_y) for coord in pair]
        self.canvas_signal.delete("all")
        self.canvas_signal.create_line(resp_points, fill="black", smooth=True, width=2)
        
        self.update_gating_calculation()

    def update_gating_calculation(self, event=None):
        """Calculate gating parameters based on respiratory signal and window settings"""
        try:
            resp_period = float(self.var_resp_period.get())
            gating_window = self.var_gating_window.get()
            
            # Generate respiratory signal for gating calculation
            resp_x = list(range(0, 360, 1))
            resp_y = self.resp_signal(resp_x)
            
            # Calculate gating threshold based on window percentage
            max_val = max(resp_y)
            threshold = max_val * (gating_window / 100.0)
            
            # Calculate T_on (time when signal is below threshold for gating)
            gating_phases = [1 if y <= threshold else 0 for y in resp_y]
            on_count = sum(gating_phases)
            
            # Convert phase count to time values
            T_on = (on_count / 360.0) * resp_period
            T_off = resp_period - T_on
            efficiency = T_on / resp_period
            
            # Store calculated parameters
            self.gating_params = {
                'T_on': T_on,
                'T_off': T_off, 
                'T_total': resp_period,
                'efficiency': efficiency
            }
            
            # Update UI labels
            self.label_ton.config(text=f"T_on: {T_on:.2f} sec")
            self.label_toff.config(text=f"T_off: {T_off:.2f} sec") 
            self.label_efficiency.config(text=f"Efficiency: {efficiency:.2f}")
            
        except ValueError:
            pass

    def f_bot_computation_with_dr(self, beam, DR=None):
        """Calculate beam-on-time with dynamic range considerations"""
        machine_db = get_current("MachineDB")
        beam_set = get_current('BeamSet')
        name = beam_set.MachineReference.MachineName
        machine = machine_db.GetTreatmentMachine(machineName=name, lockMode=None)
        bqs = machine.BeamQualities
        bq = bqs[0]
        
        # Machine parameters
        max_speed = 2000.0  # Maximum scanning speed (mm/sec)
        DR_max_limit = bq.ScanningProperties.DoseRateMaximumLimit
        DR_min_limit = 1.4  # Minimum dose rate limit
        time_resolution = 0.0001  # Time resolution for calculations

        tot_beam_MU = beam.BeamMU
        scan_times_layer = []

        # Process each segment (energy layer)
        for segment in beam.Segments:
            layer_linespacing = segment.Spots.SpotGridSpacing
            layer_weights = np.array(segment.Spots.Weights[1:])
            layer_E = segment.NominalEnergy
            Energy_ind = int(round((layer_E - 70) / 0.4))
            
            # Calculate MU per distance
            mu_per_dist = layer_weights / layer_linespacing
            
            # Calculate dose rates
            dose_rates = max_speed * mu_per_dist
            
            # Determine layer dose rate with energy-specific limits
            if len(dose_rates) > 0:
                min_dr = np.min(dose_rates)
                if min_dr > DR_max_limit[Energy_ind]:
                    layer_doserate = DR_max_limit[Energy_ind]
                elif min_dr < DR_min_limit:
                    layer_doserate = DR_min_limit
                else:
                    layer_doserate = min_dr
            else:
                layer_doserate = DR_min_limit
            
            # Apply dynamic range parameter if specified
            if DR is not None:            
                speeds = layer_doserate / mu_per_dist
                if len(speeds) > 0:
                    DR_eff = min(DR, np.max(speeds)/(1.2*np.min(speeds)))
                else:
                    DR_eff = DR
                
                max_dose_rate = np.max(dose_rates)                
                layer_doserate_internal = max_dose_rate / DR_eff
                
                # Calculate beam-on-time for each segment
                bot_internals = layer_weights / layer_doserate_internal
                speed_segments = layer_linespacing / bot_internals
                speed_mask = speed_segments > max_speed
            
                # Adjust weights for segments exceeding speed limit
                layer_weights[speed_mask] = layer_doserate_internal * layer_linespacing / max_speed
                layer_doserate = layer_doserate_internal                    

            # Calculate raw scan time
            raw_scan_time = layer_weights / layer_doserate
            
            # Round scan time to time resolution
            rounded_scan_time = time_resolution * np.round(raw_scan_time / time_resolution)
            
            scan_times_layer.append(rounded_scan_time)

        return scan_times_layer

    def calculate_total_treatment_time(self, scan_times_all_layers):
        """Calculate total treatment time considering gating parameters"""
        # Non-gated case
        if not self.gating_params or self.gating_params['T_on'] == 0:
            total_scan_time = sum([sum(layer_times) for layer_times in scan_times_all_layers])
            return {
                'total_time': total_scan_time,
                'beam_on_time': total_scan_time,
                'efficiency': 1.0
            }

        # Gated treatment calculation
        T_on = self.gating_params['T_on']
        T_total = self.gating_params['T_total']
        layer_switching_time = float(self.var_layer_switch_time.get())
        
        t = 0.0  # Total elapsed time
        total_beam_on_time = 0.0

        # Process each layer's scan times
        for i, layer_scan_times in enumerate(scan_times_all_layers):
            for scan_time in layer_scan_times:
                remaining = float(scan_time)
                
                # Distribute remaining scan time across gating cycles
                while round(remaining, 5) > 0:
                    current_cycle_position = round(t % T_total, 5)
                    if T_total - current_cycle_position < 1e-4:
                        current_cycle_position = 0.0
                    
                    # Check if currently in beam-on phase
                    if round(T_on - current_cycle_position, 5) > 0:
                        available_time = T_on - current_cycle_position
                        used_time = min(available_time, remaining)
                        remaining -= used_time
                        total_beam_on_time += used_time
                        t += used_time
                    else:
                        # Skip to next beam-on phase
                        next_on_time = t + (T_total - current_cycle_position)
                        t = next_on_time
            
            # Add layer switching time between layers
            if i < len(scan_times_all_layers) - 1:
                t += layer_switching_time

        return {
            'total_time': t,
            'beam_on_time': total_beam_on_time,
            'efficiency': total_beam_on_time / t if t > 0 else 0
        }

    def analyze_dr_range(self):
        """Perform dynamic range analysis across specified range"""
        try:
            self.status_text.config(text="Starting DR analysis...")
            self.master.update()
            
            case = get_current("Case")
            plan_name = self.var_plan_name.get()
            beam_set = case.TreatmentPlans[plan_name].BeamSets[0]

            TR = float(self.var_resp_period.get())
            
            # Get DR analysis parameters
            dr_start = int(self.var_dr_start.get())
            dr_end = int(self.var_dr_end.get())
            dr_step = int(self.var_dr_step.get())
            
            dr_values = list(range(dr_start, dr_end + 1, dr_step))
            results = {'DR': [], 'total_time': [], 'beam_on_time': [], 'efficiency': []}
            
            # Analyze each DR value
            for i, dr_val in enumerate(dr_values):
                self.status_text.config(text=f"Analyzing DR={dr_val} ({i+1}/{len(dr_values)})")
                self.master.update()
                
                # Calculate scan times for all beams at current DR
                all_layer_scan_times = []
                for beam in beam_set.Beams:
                    scan_times = self.f_bot_computation_with_dr(beam, DR=dr_val)
                    all_layer_scan_times.extend(scan_times)
                
                # Calculate total treatment time
                time_results = self.calculate_total_treatment_time(all_layer_scan_times)
                
                # Store results
                results['DR'].append(dr_val)
                results['total_time'].append(time_results['total_time'])
                results['beam_on_time'].append(time_results['beam_on_time'])
                results['efficiency'].append(time_results['efficiency'])
            
            self.dr_results = results
            self.plot_results()
            self.status_text.config(text="DR analysis completed.")
            
        except Exception as e:
            messagebox.showerror("Error", f"Analysis failed: {str(e)}")
            self.status_text.config(text="Analysis failed.")

    def data_to_canvas_coords(self, data_x, data_y, canvas_width, canvas_height, margin=50):
        """Convert data coordinates to canvas coordinates for plotting"""
        if not data_x or not data_y:
            return [], []
            
        # Calculate data ranges
        min_x, max_x = min(data_x), max(data_x)
        min_y, max_y = min(data_y), max(data_y)
        
        # Handle zero range cases
        if max_x == min_x:
            max_x = min_x + 1
        if max_y == min_y:
            max_y = min_y + 1
            
        # Convert to canvas coordinates
        graph_width = canvas_width - 2 * margin
        graph_height = canvas_height - 2 * margin
        
        canvas_x = [margin + (x - min_x) / (max_x - min_x) * graph_width for x in data_x]
        canvas_y = [canvas_height - margin - (y - min_y) / (max_y - min_y) * graph_height for y in data_y]
        
        return canvas_x, canvas_y, (min_x, max_x, min_y, max_y)

    def draw_axes_and_grid(self, canvas, canvas_width, canvas_height, data_ranges, margin=50):
        """Draw axes and grid lines on canvas"""
        min_x, max_x, min_y, max_y = data_ranges
        
        # Draw main axes
        canvas.create_line(margin, canvas_height - margin, 
                          canvas_width - margin, canvas_height - margin, 
                          fill="black", width=2)  # X-axis
        canvas.create_line(margin, margin, 
                          margin, canvas_height - margin, 
                          fill="black", width=2)  # Y-axis
        
        # Draw grid lines
        grid_lines = 9
        for i in range(grid_lines + 1):
            # Vertical grid lines
            x = margin + i * (canvas_width - 2 * margin) / grid_lines
            canvas.create_line(x, margin, x, canvas_height - margin, 
                              fill="lightgray", width=1)
            
            # Horizontal grid lines
            y = margin + i * (canvas_height - 2 * margin) / grid_lines
            canvas.create_line(margin, y, canvas_width - margin, y, 
                              fill="lightgray", width=1)
        
        # Draw axis labels
        for i in range(grid_lines + 1):
            # X-axis labels
            x_val = min_x + i * (max_x - min_x) / grid_lines
            x_pos = margin + i * (canvas_width - 2 * margin) / grid_lines
            canvas.create_text(x_pos, canvas_height - margin + 20, 
                              text=f"{x_val:.0f}", anchor="n")
            
            # Y-axis labels
            y_val = min_y + i * (max_y - min_y) / grid_lines
            y_pos = canvas_height - margin - i * (canvas_height - 2 * margin) / grid_lines
            canvas.create_text(margin - 10, y_pos, 
                              text=f"{y_val:.0f}", anchor="e")

    def draw_data_points(self, canvas, canvas_x, canvas_y):
        """Draw data points and connecting lines on canvas"""
        if len(canvas_x) < 2:
            return
            
        # Draw connecting lines
        line_points = [coord for pair in zip(canvas_x, canvas_y) for coord in pair]
        canvas.create_line(line_points, fill="blue", width=2, smooth=True)
        
        # Draw data point markers
        for x, y in zip(canvas_x, canvas_y):
            canvas.create_oval(x-3, y-3, x+3, y+3, fill="red", outline="darkred")

    def draw_labels(self, canvas, canvas_width, canvas_height):
        """Draw axis titles and graph title"""
        # X-axis title
        canvas.create_text(canvas_width // 2, canvas_height - 10, 
                          text="Dynamic Range", anchor="s", font=("Arial", 10, "bold"))
        
        # Y-axis title
        canvas.create_text(15, canvas_height // 2, 
                          text="Tx Time (sec)", anchor="center", 
                          font=("Arial", 10, "bold"))
        
        # Graph title
        canvas.create_text(canvas_width // 2, 20, 
                          text="DR vs Treatment Time for every fields", anchor="n", 
                          font=("Arial", 12, "bold"))

    def plot_results(self):
        """Plot results using canvas graphics"""
        if not self.dr_results or not self.dr_results['DR']:
            self.canvas_graph.delete("all")
            return
        
        # Clear canvas and schedule delayed drawing
        self.canvas_graph.delete("all")
        self.canvas_graph.update_idletasks()
        self.canvas_graph.update()
        
        self.master.after(10, self._draw_graph_delayed)

    def _draw_graph_delayed(self):
        """Delayed graph drawing to ensure proper canvas sizing"""
        # Get canvas dimensions
        canvas_width = self.canvas_graph.winfo_width()
        canvas_height = self.canvas_graph.winfo_height()
        
        # Ensure minimum size
        if canvas_width < 100:
            canvas_width = 800
        if canvas_height < 100:
            canvas_height = 400
        
        # Draw complete graph
        canvas_x, canvas_y, data_ranges = self.data_to_canvas_coords(
            self.dr_results['DR'], 
            self.dr_results['total_time'],
            canvas_width, 
            canvas_height
        )
        
        self.draw_axes_and_grid(self.canvas_graph, canvas_width, canvas_height, data_ranges)
        self.draw_data_points(self.canvas_graph, canvas_x, canvas_y)
        self.draw_labels(self.canvas_graph, canvas_width, canvas_height)
        
        self.canvas_graph.update()

    def export_results(self):
        """Export analysis results to CSV file"""
        if not self.dr_results:
            messagebox.showwarning("Warning", "No results to export.")
            return
            
        filename = filedialog.asksaveasfilename(defaultextension=".csv",
                                               filetypes=[("CSV files", "*.csv")])
        if filename:
            import csv
            with open(filename, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['DR', 'Total_Time', 'Beam_On_Time', 'Efficiency'])
                for i in range(len(self.dr_results['DR'])):
                    writer.writerow([self.dr_results['DR'][i], 
                                   self.dr_results['total_time'][i],
                                   self.dr_results['beam_on_time'][i],
                                   self.dr_results['efficiency'][i]])
            messagebox.showinfo("Success", f"Results exported to {filename}")

if __name__ == "__main__":
    patient = get_current("Patient")
    case = get_current('Case')

    root = tk.Tk()
    app = GUIProgram(root)
    root.mainloop()