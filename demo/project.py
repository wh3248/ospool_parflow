"""
Create and populate a parflow directory with input files and a runscript to be
used by parflow to generate output simulation files for a subset domain.
"""

# pylint: disable = C0301,R0914,W0632
import os
import shutil
import datetime
import parflow
import numpy as np
import hf_hydrodata as hf
import subsettools as st


def create_project(project_options: dict, directory_path: str = "project_dir") -> str:
    """
    Create a parflow project directory to collect input files about a domain
    that may be used to run parflow to generate the output files of a simulation.

    Parameters:
        directory_path:     Path name to a directory where the parflow files are created.
        project_options:    A dict of keys with options to run a parflow simulation.

    The project_options dict supports the keys:
        run_type:       Either "transient" or "spinup" (defaults to "transient").
        template:       The path to a parflow yaml file for parflow parameters (optional).
        start_date:     The start date of the run as as string YYYY-mm-dd.
        end_date:       The end date of the run as as string YYYY-mm-dd.
        time_steps:     The number of timesteps to run parflow (defaults to hours between start/end).
        huc_id:         An array or comma seperated list of HUC id to subset inputs and outputs (optional).
        grid_bounds:    The conus2 points to subset (min_x, min_y, max_x, max_y) (optional).
        latlon_bounds:  The latlon bounds to subset ((min,lat, min_lon),(max_lat, max_lon) (optional).
        forcing_day:    Use fixed forcing data for every input hour using this day (YYYY-mm-dd).
        forcing_precip: Use this fixed precipitation value for every input hour (optional).
        grid:           The grid size (only conus2 is supported now) (defaults to conus2).
        topology:       An array or tuple [p, q, r] that defines the topology for generated pfb files.

    Only one of hucs, grid_bounds or latlon_bounds may be provided.
    If template is provided this overrides the run_type.
    The topology defaults to (1,1,1).

    The hucs may be a string of a comma seperated list of HUC id or an array of HUC id.

    Collects all required parflow input files into the directory_path.
    This uses subsettools and hf_hydrodata to subset the input files to the domain
    defined by hucs, grid_bounds, or latlon_bounds.

    This uses pf_tools to run a parflow simulation and writes the output simulated files
    into the same directory_path.

    For more advanced usage of parflow you may specify your own parflow template and collect
    additional input files into the project directory used by your template.

    Returns:
        The runscript as a the path to the parflow yaml file generated in the directory path.

    Example:

    .. code-block:: python
        project_options = {
            "run_type": "transient",
            "grid_bounds": [3749, 1583, 3759, 1593],
            "start_date": "2005-10-01",
            "end_date": "2005-10-02",
            "time_steps": 10,
            "topology": (1, 1, 1)
        }
        directory_path = "./trivial"

        # Create the parflow model and generated input files
        runscript_path = project.create_project(project_options, directory_path)

        # Run the parflow model
        model = parflow.Run.from_definition(runscript_path)
        model.run()

    """

    runname = os.path.basename(directory_path)
    run_type = project_options.get("run_type")
    template = project_options.get("template")
    if not template:
        if run_type == "transient":
            template = "conus2_transient_solid.yaml"
        elif run_type == "spinup":
            template = "conus2_spinup_solid.yaml"
        elif run_type:
            raise ValueError(
                f"Unsupported run_type '{run_type}'. Must be transient or spinup."
            )
        else:
            template = "conus2_transient_solid.yaml"
    else:
        template = "conus2_transient_solid.yaml"

    runscript_path = _create_runscript(runname, directory_path, template)
    _create_topology(runscript_path, project_options)
    _create_static_and_forcing(runscript_path, project_options, runname)
    _create_dist_files(runscript_path, project_options)

    return runscript_path


def _create_runscript(
    runname: str,
    directory_path: str,
    template_path="conus2_transient_solid.yaml",
):
    """
    Create a parflow model using the template.
    Returns:
        the path to the runscript of the model
    """

    directory_path = os.path.abspath(directory_path)

    os.makedirs(directory_path, exist_ok=True)
    template_dir = os.path.dirname(os.path.abspath(__file__))
    if not template_path.startswith("/"):
        template_path = os.path.join(template_dir, template_path)
    parflow.tools.settings.set_working_directory(directory_path)
    runscript_path = os.path.abspath(f"{directory_path}/{runname}.yaml")

    # Create Parfow runscript_path using the template if it does not exist yet
    if not os.path.exists(runscript_path):
        shutil.copy(template_path, runscript_path)
        model = parflow.Run.from_definition(runscript_path)
        model.write(file_format="yaml")

    return runscript_path


def _create_topology(runscript_path: str, project_options: dict):
    """
    Create the topology files and add the references to the model and runscript.yaml file
    """
    model = parflow.Run.from_definition(runscript_path)
    _, grid, ij_bounds, latlon_bounds, _, _ = _get_time_space_options(project_options)

    topology = project_options.get("topology")
    topology = list(topology) if isinstance(topology, tuple) else topology
    topology = [1, 1, 1] if not topology else topology
    if not len(topology) == 3:
        raise ValueError(
            "The topology option in project options must be an array [p, q, r]"
        )

    p = int(topology[0])
    q = int(topology[1])
    r = int(topology[2])
    if r != 1:
        raise ValueError("The r dimension of the topology must be 1")
    model.Process.Topology.P = p
    model.Process.Topology.Q = q
    model.Process.Topology.R = r
    model.FileVersion = 4

    ij_bounds, _ = st.define_latlon_domain(latlon_bounds, grid)

    model.ComputationalGrid.Lower.X = ij_bounds[0]
    model.ComputationalGrid.Lower.Y = ij_bounds[1]
    model.ComputationalGrid.Lower.Z = 0.0

    # Define the size of each grid cell. The length units are the same as those on hydraulic conductivity, here that is meters.
    model.ComputationalGrid.DX = 1000.0
    model.ComputationalGrid.DY = 1000.0
    model.ComputationalGrid.DZ = 200.0

    # Define the number of grid blocks in the domain.
    model.ComputationalGrid.NX = ij_bounds[2] - ij_bounds[0]
    model.ComputationalGrid.NY = ij_bounds[3] - ij_bounds[1]
    if grid == "conus1":
        model.ComputationalGrid.NZ = 5
    elif grid == "conus2":
        model.ComputationalGrid.NZ = 10

    model.write(file_format="yaml")


def _create_static_and_forcing(runscript_path: str, options: dict, runname: str):
    """
    Create the static input and forcing files and add the references to the model and runscript.yaml file
    """
    model = parflow.Run.from_definition(runscript_path)
    directory_path = os.path.dirname(runscript_path)

    mask, grid, ij_bounds, _, start_date, end_date = _get_time_space_options(options)

    st.write_mask_solid(mask=mask, grid=grid, write_dir=directory_path)

    var_ds = "conus2_domain"
    static_paths = st.subset_static(ij_bounds, dataset=var_ds, write_dir=directory_path)
    st.config_clm(
        ij_bounds,
        start=start_date,
        end=end_date,
        dataset=var_ds,
        write_dir=directory_path,
    )

    forcing_dir_path = directory_path
    os.makedirs(forcing_dir_path, exist_ok=True)
    forcing_day = options.get("forcing_day", None)
    forcing_ds = "CW3E"
    if forcing_day:
        # use fixed values for all forcing hour inputs
        precip = options.get("precip", None)
        start_time_dt = datetime.datetime.strptime(start_date, "%Y-%m-%d")
        end_time_dt = datetime.datetime.strptime(end_date, "%Y-%m-%d")
        forcing_variables = [
            "downward_shortwave",
            "precipitation",
            "downward_longwave",
            "specific_humidity",
            "air_temp",
            "atmospheric_pressure",
            "east_windspeed",
            "north_windspeed",
        ]
        for variable in forcing_variables:
            # Get the forcing data for all variables
            options = {
                "dataset": forcing_ds,
                "variable": variable,
                "grid_bounds": list(ij_bounds),
                "temporal_resolution": "daily",
                "start_time": forcing_day,
                "aggregation": "sum" if variable == "precipitation" else "mean",
                "dataset_version": "1.0",
            }
            metadata = hf.get_catalog_entry(options)
            dataset_var = (
                "Press"
                if variable == "atmospheric_pressure"
                else "Temp" if variable == "air_temp" else metadata.get("dataset_var")
            )

            data = hf.get_gridded_data(options)
            if variable == "precipitation":
                data = data / 24
                if precip:
                    data[:, :, :] = float(precip)
            day_data = np.zeros((24, data.shape[1], data.shape[2]))
            for i in range(0, 24):
                # Set the data to be the same for all 24 hours in the PFB file
                day_data[i, :, :] = data[0, :, :]
            dt = start_time_dt
            day = 1
            while dt < end_time_dt:
                # Create hourly pfb files for each day in the parflow run range to all be the same
                forcing_file_path = f"{forcing_dir_path}/{forcing_ds}.{dataset_var}.{day:06d}_to_{day+23:06d}.pfb"
                parflow.write_pfb(forcing_file_path, day_data)
                dt = dt + datetime.timedelta(days=1)
                day = day + 24
    else:
        # Get the forcing data from the CW3E dataset for the days in the parflow run range
        st.subset_forcing(
            ij_bounds,
            grid=grid,
            start=start_date,
            end=end_date,
            dataset=forcing_ds,
            write_dir=forcing_dir_path,
        )

    # Update the runscript yaml file with the forcing_dir_path
    st.edit_runscript_for_subset(
        ij_bounds,
        runscript_path=runscript_path,
        runname=runname,
        forcing_dir=forcing_dir_path,
    )

    # Update the file names of the generated parflow static files
    init_press_path = os.path.basename(static_paths["ss_pressure_head"])
    depth_to_bedrock_path = os.path.basename(static_paths["pf_flowbarrier"])

    st.change_filename_values(
        runscript_path=runscript_path,
        init_press=init_press_path,
        depth_to_bedrock=depth_to_bedrock_path,
    )

    # Set the forcing file dataset name in the runscript yaml file
    model = parflow.Run.from_definition(runscript_path)
    model.Solver.CLM.MetFileName = "CW3E"
    model.write(file_format="yaml")


def _create_dist_files(runscript_path: str, options: dict):
    """
    Create the parflow .dist files for the generated pfb files in the parflow directory.
    """
    p = int(options.get("p", "1"))
    q = int(options.get("q", "1"))

    st.dist_run(
        topo_p=p,
        topo_q=q,
        runscript_path=runscript_path,
        dist_clim_forcing=True,
    )
    model = parflow.Run.from_definition(runscript_path)

    # Set the timesteps to use in the parflow run
    time_steps = options.get("time_steps", None)
    if time_steps is None:
        # If time_steps is not set in the options use the hours between start and end time
        _, _, _, _, start_date, end_date = _get_time_space_options(options)
        start_time_dt = datetime.datetime.strptime(start_date, "%Y-%m-%d")
        end_time_dt = datetime.datetime.strptime(end_date, "%Y-%m-%d")
        days_between = (end_time_dt - start_time_dt).days
        model.TimingInfo.StopTime = 24 * int(days_between)
    else:
        # If time_steps is set in the options then use that number of steps
        model.TimingInfo.StopTime = int(time_steps)

    # Reset the NZ that can be incorrectly set by st.dist_run
    model.ComputationalGrid.NZ = 10
    model.write(file_format="yaml")


def _get_time_space_options(options):
    """
    Get the time and space options from the input options.
    Returns:
        (grid, ij_bounds, latlon_bounds, start_date, end_date)
    """

    grid_bounds = options.get("grid_bounds", None)
    latlon_bounds = options.get("latlon_bounds", None)
    huc_id = options.get("huc_id", None)
    grid = options.get("grid", "conus2")
    start_date = options.get("start_date", "2001-01-01")
    end_date = options.get("end_date", "2001-01-02")
    if huc_id:
        hucs = list(huc_id) if isinstance(huc_id, tuple) else huc_id if isinstance(huc_id, list) else huc_id.split(",")
        ij_bounds, mask = st.define_huc_domain(hucs=hucs, grid=grid)
        lat_min, lon_min = hf.to_latlon(grid, ij_bounds[0], ij_bounds[1])
        lat_max, lon_max = hf.to_latlon(grid, ij_bounds[2] - 1, ij_bounds[3] - 1)
        latlon_bounds = [[lat_min, lon_min], [lat_max, lon_max]]
    elif grid_bounds:
        lat_min, lon_min = hf.to_latlon(grid, grid_bounds[0], grid_bounds[1])
        lat_max, lon_max = hf.to_latlon(grid, grid_bounds[2] - 1, grid_bounds[3] - 1)
        latlon_bounds = [[lat_min, lon_min], [lat_max, lon_max]]
        ij_bounds, mask = st.define_latlon_domain(latlon_bounds, grid)
    elif latlon_bounds:
        if len(latlon_bounds) != 2:
            raise ValueError("The latlon_bounds must be an array of 2 lat/lon pairs")
        if len(latlon_bounds[0]) != 2:
            raise ValueError("The latlon_bounds must be an array of 2 lat/lon pairs")
        ij_bounds, mask = st.define_latlon_domain(latlon_bounds, grid)
    else:
        raise ValueError("Must specify in options hucs, grid_bounds, or latlon_bounds")
    return (mask, grid, ij_bounds, latlon_bounds, start_date, end_date)
