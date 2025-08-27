from shiny import App, reactive, render, ui, module
from shinywidgets import output_widget, render_widget
import tempfile
import os

@module.ui
def query_heatmap_ui():
    return (
        ui.output_text("name_text"),
        output_widget("generated_plot"),
        ui.row(
            ui.download_link("download_data", "Download PNG"),  # Add download button
            ui.input_action_link("delete", "Delete")
        )
    )

@module.server
def query_heatmap_server(input, output, session, remove_id, plot, name: str):
    @render_widget
    def generated_plot():
        return plot
    
    @render.text
    def name_text():
        return name

    @reactive.effect
    @reactive.event(input.delete)
    def _():
        ui.remove_ui(selector=f"#{remove_id}")

    @render.download(filename="heatmap.png", media_type="image/png")  # Define the download handler
    def download_data():
        f = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
        print("[DEBUG] Saving heatmap to:", f.name)
        try:
            print(f"Plot type: {type(plot)}")
            plot.write_image(f.name)  # Use kaleido to write the PNG
            f.flush()
            file_size = os.path.getsize(f.name)
            print(f"[DEBUG] PNG written: {f.name}, size: {file_size} bytes")
        finally:
            f.close()

        yield f.name  # Yield the file path to trigger download
