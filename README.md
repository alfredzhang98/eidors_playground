# EIDORS Interactive Simulation Interface

A MATLAB GUI for interactive electrical impedance tomography (EIT) simulation and reconstruction using the EIDORS toolkit.  Users can place circular inclusions in a 2D tank model, adjust inclusion properties, and perform both standard and “high-quality” reconstructions, all in an easy-to-use visual interface.

---

## Features

- **Interactive inclusion placement**  
  Click inside the tank boundary to add circular inclusions at any position.

- **Adjustable inclusion parameters**  
  - **Radius** (0.5 – 4.0 UI units)  
  - **Conductivity multiplier** (0.1 – 10.0, default 5.0)  
  - **Regularization parameter** for inverse solver (1e-4 – 1e0)

- **Real-time reconstruction**  
  Toggle “Real-time Recon” to automatically re-solve after each placement.

- **Standard & High-Quality Reconstruction**  
  - **Standard**: single solver (`inv_solve_diff_GN_one_step`) with Gaussian HPF prior  
  - **High-Quality**: tests multiple priors (HPF, TV, Laplace, NOSER), scores results, and displays the best reconstruction

- **Model visualization**  
  - Show forward (high-res) and reconstruction (coarse) mesh  
  - Display the true conductivity distribution of placed inclusions

- **Scrollable log panel**  
  A `listbox` control shows timestamped log entries and auto-scrolls to the latest message.

---

## Requirements

- MATLAB R2024a or later  
- [EIDORS v3.x](https://eidors3d.sourceforge.net/) (add its `startup_eidors` folder to your MATLAB path)  
- Java Swing (built into MATLAB)  
- `findjobj` (optional; only if you revert to an `edit` box for logs)

---

## Installation

1. **Clone or download** this repository into your MATLAB working folder.  
2. In MATLAB, **add** both this folder and the EIDORS startup script to your path:  
   ```matlab
   addpath(genpath('path/to/eidors'));
   startup_eidors;
   addpath(genpath('path/to/EIDORS-Interactive-UI'));
````

3. **Verify** EIDORS is working by running:

   ```matlab
   which eidors_obj
   ```

   You should see the path to `eidors_obj.m`.

---

## Usage

From the MATLAB command window, simply run:

```matlab
eidors_interactive_ui
```

The main window will open:

1. **Left Panel:**

   * Click inside the circular tank to place inclusions.
   * Inclusions are numbered in placement order.

2. **Control Panel (upper right):**

   * **Inclusion Radius**: adjust size before clicking.
   * **Conductivity Mult**: default `5.0`, relative to background conductivity.
   * **Regularization**: default `1e-1`, controls solver smoothness.
   * **Reconstruct**: run a one-step reconstruction.
   * **High Quality Recon**: test multiple priors and pick the best result.
   * **Real-time Recon**: automatically reconstruct after each placement.
   * **Show Forward Model / Recon Model / True Distrib**: visualize meshes and ground truth.

3. **Results Panel (lower right):**

   * Displays reconstructed conductivity map or true inclusion distribution.
   * Colorbar and dynamic `CLim` based on data range.

4. **Log Panel (bottom of Control Panel):**

   * Shows time-stamped messages.
   * Auto-scrolls to newest entry.

---

## Customization

* **Background Conductivity**
  In `init_app_data`, change:

  ```matlab
  data.homg_img = mk_image(data.fmdl_fwd, 5);
  ```

  to any desired value.

* **Default Conductivity Slider**
  In `create_controls`, set slider `Value` and `String`:

  ```matlab
  'Value', 5,  % default conductivity
  ...
  'String', '5.0'
  ```

* **Log Control Type**
  If you prefer the old multiline `edit` box, you can revert and use `findjobj` + `setCaretPosition` in `update_info`, but the `listbox` version is more robust.

---

## File Overview

* **eidors\_interactive\_ui.m**
  Main entry point—defines GUI layout, callbacks, and core logic.
* **init\_app\_data**
  Loads forward/recon models, sets up uniform background image/data.
* **create\_ui / create\_controls**
  Constructs panels, sliders, buttons, and log control.
* **update\_info**
  Appends timestamped log entries and auto-scrolls.
* **perform\_reconstruction**
  Standard reconstruction workflow with debug logging.
* **advanced\_reconstruction**
  High-quality multi-prior solver testing and selection.
* **display / plotting helpers**
  Mesh display, slice rendering, and statistics reporting.

---

## Troubleshooting

* **“Please initialize EIDORS first”**
  Ensure you’ve run `startup_eidors` and added its path.
* **Slow performance**
  High-res forward model can be slow; use “Real-time Recon” sparingly or reduce mesh density.
* **Log not auto-scrolling**
  Confirm you’re using the `listbox` version of `text_info`. Check that `ListboxTop` is set in `update_info`.

---

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.

---

Enjoy interactive EIT simulations with the power of EIDORS!
Feedback and contributions welcome.
