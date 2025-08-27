from shiny import App, reactive, render, ui, module
import db

@module.ui
def query_output_ui():
    return (
        ui.output_text("name_text"),
        ui.input_text("filter", "Filter Results", autocomplete=False),  # Add label for clarity
        ui.output_data_frame("results"),
        ui.row(
            ui.download_link("download_data", "Download CSV"),  # Add download button
            ui.input_action_link("delete", "Delete")
        )
    )

@module.server
def query_output_server(input, output, session, remove_id, result, name: str):
    safe_name = name.replace(" ", "_").replace(",", "").replace("=", "-")

    @render.data_frame
    def results():
        return results_filtered()
    
    @render.text
    def name_text():
        return name
    
    @reactive.calc
    def results_filtered():
        filter_text = input.filter().strip().lower()
        if filter_text:
            filtered_result = result[result['Gene'].str.lower().str.contains(filter_text)]
            return filtered_result
        else:
            return result

    @reactive.effect
    @reactive.event(input.delete)
    def _():
        ui.remove_ui(selector=f"#{remove_id}")

   # @render.download(filename="table.csv")  # Define the download handler
   # def download_data():
   #     yield db.df_to_csv(result)  # Convert DataFrame to CSV for download
   #     return

    @render.download(filename=lambda: f"{safe_name}.csv")  # Use name as the filename
    def download_data():
        yield db.df_to_csv(result)

     

