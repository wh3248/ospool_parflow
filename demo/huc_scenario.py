import project
import parflow
import hf_hydrodata as hf
import time

def main():
    """
    Test generating a parflow directory and execute model and assert start/end pressure values for a box
    Use a box with radius 5 around the same target point that is the center of HUC 02080203
    This should get the same answer as the HUC test, but with a smaller parflow domain.
    """

    try:
        duration_start = time.time()
        hf.register_api_pin("wh3248@princeton.edu", "0000")
        runname = "trival"
        directory_path = f"./{runname}"

        start_date = "2005-10-01"
        end_date = "2005-10-02"
        target_x = 3754
        target_y = 1588
        target_radius = 5
        time_steps = 10
        grid_bounds = [
            target_x - target_radius,
            target_y - target_radius,
            target_x + target_radius,
            target_y + target_radius,
        ]
        parflow_options = {
            "grid_bounds": grid_bounds,
            "grid": "conus2",
            "start_date": start_date,
            "end_date": end_date,
            "time_steps": time_steps,
            "forcing_day": start_date,
        }

        # Create the parflow model and generated input files
        runscript_path = project.create_project(parflow_options, directory_path, )
        model = parflow.Run.from_definition(runscript_path)
        model.write(file_format="yaml")

        input_duration = time.time() - duration_start

        parflow_start = time.time()
        # Run the parflow model
        model.run()
        parflow_duration = time.time() - parflow_start
        full_duration = time.time() - duration_start
        print(f"Duration {input_duration} seconds to collect inputs.")
        print(f"Duration {parflow_duration} seconds to run parflow.")
        print(f"Duration {full_duration} seconds to collect inputs and run parflow.")
        print()

        # Verify the result
        time_step = time_steps - 1
        nx = model.ComputationalGrid.NX
        ny = model.ComputationalGrid.NY
        nz = model.ComputationalGrid.NZ
        x = int(nx / 2)
        y = int(ny / 2)
        z = nz - 1
        out_path = f"{directory_path}/{runname}.out.press.{time_step:05d}.pfb"
        out_press_np = parflow.read_pfb(out_path)
        print(f"OUT PRESS ({z},{y},{x}) {out_press_np[z, y, x]} [{time_step}]")
        top_layer_pressure = out_press_np
        assert round(top_layer_pressure[z, y, x], 5) == 0.00325


    except Exception as e:
        raise e
    
main()