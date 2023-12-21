"""MLCube handler file"""
import typer
import subprocess


app = typer.Typer()


def exec_python(cmd: str) -> None:
    """Execute a python script as a subprocess

    Args:
        cmd (str): command to run as would be written inside the terminal
    """
    splitted_cmd = cmd.split()
    process = subprocess.Popen(splitted_cmd, cwd=".")
    process.wait()
    assert process.returncode == 0, f"command failed: {cmd}"


@app.command("prepare")
def prepare(
    data_path: str = typer.Option(..., "--data_path"),
    labels_path: str = typer.Option(..., "--labels_path"),
    parameters_file: str = typer.Option(..., "--parameters_file"),
    models_path: str = typer.Option(..., "--models"),
    output_path: str = typer.Option(..., "--output_path"),
    output_labels_path: str = typer.Option(..., "--output_labels_path"),
    report_file: str = typer.Option(..., "--report_file"),
    metadata_path: str = typer.Option(..., "--metadata_path"),
):
    cmd = f"python3 project/prepare.py --data_path={data_path} --labels_path={labels_path} --models_path={models_path} --data_out={output_path} --labels_out={output_labels_path} --report={report_file} --parameters={parameters_file} --metadata_path={metadata_path}"
    exec_python(cmd)


@app.command("sanity_check")
def sanity_check(
    data_path: str = typer.Option(..., "--data_path"),
    labels_path: str = typer.Option(..., "--labels_path"),
    parameters_file: str = typer.Option(..., "--parameters_file"),
    metadata_path: str = typer.Option(..., "--metadata_path"),
):
    # Modify the sanity_check command as needed
    cmd = f"python3 project/sanity_check.py --data_path={data_path} --labels_path={labels_path} --metadata={metadata_path}"
    exec_python(cmd)


@app.command("statistics")
def sanity_check(
    data_path: str = typer.Option(..., "--data_path"),
    labels_path: str = typer.Option(..., "--labels_path"),
    parameters_file: str = typer.Option(..., "--parameters_file"),
    metadata_path: str = typer.Option(..., "--metadata_path"),
    out_path: str = typer.Option(..., "--output_path"),
):
    # Modify the statistics command as needed
    cmd = f"python3 project/statistics.py --data_path={data_path} --labels_path={labels_path} --out_file={out_path} --metadata={metadata_path}"
    exec_python(cmd)


if __name__ == "__main__":
    app()
