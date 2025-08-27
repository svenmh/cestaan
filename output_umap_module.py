from shiny import module, ui, reactive, render
import os

@module.ui
def query_umap_ui():
    return ui.div(

        ui.output_text("name_text"),
        ui.output_ui("output_umap_test"),
        ui.row(
            #ui.download_link("download_data", "Download PNG"),  # Add download button
            ui.column(2),
            ui.column(4, ui.output_ui("download_href")),
            ui.column(4, ui.input_action_link("delete", "Delete"))
        )
        
     )
    

@module.server
def query_umap_server(input, output, session, remove_id, img_path, name: str, is_two_genes=False):
    print(f"[DEBUG 2A] Entered query_umap_server with id={remove_id}, path={img_path}, output={output}")


    @render.ui
    def output_umap_test():
        print("[INFO 2A] entered output_umap function")
        full_path = os.path.abspath(os.path.join("www", img_path))
        file_exists = os.path.exists(full_path)
        print(f"[DEBUG] full_path: {full_path}, exists: {file_exists}")

        rel_path = os.path.basename(img_path)
        print(f"[query_umap_server] Rendering image from: {img_path} -> /{rel_path}")
        print("[DEBUG 2B] Reached generated_plot render function.")

        img_style = "max-width: 60%; border: 1px solid #ccc;"
        if is_two_genes:
            img_style = "max-width: 90%; border: 1px solid #ccc;"

        return ui.tags.div(
                #ui.tags.p(f"Image source: /{img_path}"),
                #ui.tags.p(f"File exists: {file_exists}"),
                ui.tags.img(src=f"/{img_path}"
                    , style= img_style)
                    )
        #return ui.tags.img(src=f"/{rel_path}", style="max-width: 100%; border: 1px solid #ccc;")
    #print("[DEBUG 2BA] Forced evaulation of output_umap:", output_umap())

    @render.text
    def name_text():
        print("[INFO 2C] entered name_text function")
        return name
    
    @render.ui
    def download_href():
        download_filename = os.path.basename(img_path)
        return ui.tags.a(
                "Download PNG",
                href=f"/{img_path}",
                download=download_filename,
                target="_blank",
                style="margin-right: 5px;"
            )

    @reactive.effect
    @reactive.event(input.delete)
    def _():
        print(f"[INFO 2D] Removing UI for {remove_id}")
        ui.remove_ui(selector=f"#{remove_id}")

        # Required to activate render functions
    
    print("[DEBUG 2E] About to assign outputs.")
    output.output_umap_test = output_umap_test
    output.name_text = name_text
    output.download_href = download_href

    print("[DEBUG 2F] Render functions assigned.")

