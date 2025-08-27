import faicons as fa #fontawesome icons for app usage
import plotly.express as px
from shinywidgets import output_widget, render_plotly
from shiny import App, reactive, render, ui, module
import db
from output_heatmap_module import query_heatmap_server, query_heatmap_ui
from output_table_module import query_output_server, query_output_ui
from output_umap_module import query_umap_server, query_umap_ui
import heatmap
from pathlib import Path
import subprocess
import os
from pdf2image import convert_from_path
from PIL import Image
import uuid
from datetime import datetime
import time

CESTAAN_ROOT=os.environ['PWD']

TEST_IMAGE_PATH = "umap_d94499a4129f4fc580709eb9aa88cacd.png"
here = Path(__file__).parent #defines our current dicrectory
current_umap_path = reactive.Value(None)
static_dir = Path(__file__).parent / "www"


#UI definitions for the site
app_ui = ui.page_sidebar(
    
    #sidebar components
    ui.sidebar(
        #Lab logo on top of sidebar as clickable link to lab website
        ui.tags.a(
            ui.tags.img(
                src="lab_logo.png",
                style="display: block; margin-left: auto; margin-right: auto; width: 100%; max-width: 180px; opacity: 0.9;"
                ),
            href="https://murphylab.princeton.edu/",
            target="_blank",
            style="display:block;"
            ),
        ui.h3("Options"), #Header
        ui.input_select("genotype", "Genotype", {"n2": "N2", "daf2": "daf-2","Male":"Male", "Herm": "him-8 Herms"}), #Genotype dropdown menu
        ui.div(style="flex-grow: 1;"),  # Spacer that fills the entire height
        
        ui.div(
            ui.input_dark_mode(mode="light"), 
            style="text-align: right; margin-bottom: 1rem;"
            )

    ),

    #Add CSS styling
    ui.include_css("style.css"),

    #Main title of the site --  outside of the navigation pane
    ui.h2(
        "CeSTAAN: ",
        ui.span("C. elegans", style ="font-style: italic;"),
        " Single nucleus Transcriptomic Atlas of Adult Neurons"
        ),

    #Main content area with tabs
    ui.page_fillable(
        ui.navset_tab(   #Tab navigation bar
            ui.nav_panel("Home",
                ui.page_fluid(
                    ui.br(), #Adds a verticle space between tabs and title
                    ui.h3("Welcome to Our Adult Single Nucleus RNA Sequencing Data"), #Header
                    ui.layout_columns(
                        ui.card(
                            ui.div(
                                ui.div(
                                    ui.output_image("umap_image"),
                                    style="""padding: 4rem;text-align: center; background-color: rgba(255,255,255,0.85); border-radius: 2rem;vertical-align:text-top; width: 100%; max-width: 100%;"""
                                    ),
                                style="text-align: center; margin-bottom: 2rem;"
                            )
                        ),
                        ui.card(
                            ui.h4("Overview of This Site"),
                            ui.tags.p(
                                "The data on this website and behavioral analysis and findings come from our publications:"
                                ),
                            ui.tags.a(
                                 "St. Ange and Weng et al. Cell Genomics 2024.",
                                 href="https://doi.org/10.1016/j.xgen.2024.100720",
                                 target="_blank",
                                 style="display:block; margin-bottom:0rem; font-weight:bold; font-size:1.1rem; color:#0066cc;"
                                ),
                            ui.tags.a(
                                "Morillo et al. Cell Reports 2025.",
                                href="https://doi.org/10.1016/j.celrep.2025.116016",
                                target="_blank",
                                style="display:block;margin-top:0rem; margin-bottom:1rem; font-weight:bold; font-size:1.1rem; color:#0066cc;"
                                                                                                                                                                                                    ),
                            ui.tags.ul(
                                ui.tags.li(ui.tags.b("What is our dataset? "), "We developed a protocol to isolate and sequence individual neuronal nuclei from adult ", ui.tags.i("C. elegans"), ". Here we provide an easy way to access our high resolution transcriptomic data for adult animals' neurons."),
                                ui.tags.li(ui.tags.b("Why single nucleus instead of single cell? "), "In certain tissue types, the cell isolation protocol causes shearing and breakage of the desired cells. In neurons, this includes axon breakage and cell leakage. This makes it difficult to capture synaptically localized mRNA transcripts for example."),
                                ui.tags.li(ui.tags.b("What genotypes have you sequenced? "), "We have sequenced Adult Day 1 N2 and ", ui.tags.i("daf-2 "), "animals as well as adult males and genotypically and age matched hermaphrodites.")
                                ),


                   # ui.div(
                       # ui.div(
                           # ui.output_image("umap_image"), #Placeholder to show imaged
                           # style="""
                               # width: 100%;
                               # max-width: 100%;
                               # background-color:rgba(255, 255, 255, 0.85);
                               # padding: 5rem;
                               # border-radius: 1rem;
                               # text-align: center;
                               # justify-content: center;
                               # display: inline-flex;
                           # """
                       # ),
                       # style="text-align: center; margin-bottom: 4rem;"
                            
                        style="""
                        padding:2rem;
                        font-size: 1.1rem;
                        line-height:1.6;
                        """ 
                        ),
                    col_widths =[6,6]
                    ),
                style ="margin-bottom:4rem;"
                )
            ),
            # --- Tab: Gene Expression by celltype ---
            ui.nav_panel("Gene expression by cell type",
                ui.row(
                    ui.column(
                        3,
                        ui.br(),
                        ui.h3("Query"),  # Updated to a static label
                        ui.input_selectize("cell_type", "Select cell type", choices=[]),
                        ui.input_select("threshold", "Select threshold", ["1", "2", "3", "4"]),
                        ui.input_action_button("add_query_gene_expr", "Run Query"),
                    )
                )
            ),

            # --- Tab: Cell type by gene expression ---
            ui.nav_panel("Cell type by gene expression",
                ui.row(
                    ui.column(
                        3,
                        ui.br(),
                        ui.h3("Query"),
                        ui.input_text("gene_name", "Select gene name", placeholder="Ex. daf-7"),
                        ui.input_select("gene_threshold", "Select threshold", ["1", "2", "3", "4"]),
                        ui.input_action_button("add_query_cell_type", "Run Query"),
                    ),
                    ui.column(2),

                    # ---Bulk Download Area ---
                    ui.column(
                        3,
                        ui.br(),
                        ui.h3("Bulk Download"),  # Updated to a static label
                        ui.input_text_area("gene_bulk_names", "Select all genes to download", placeholder="List each gene on its own line. You can also use an asterisk to get all genes with a prefix (ex. flp-*)", spellcheck="false", height="100px"),
                        ui.input_select("gene_bulk_threshold", "Select threshold", ["1", "2", "3", "4"]),
                        ui.download_button("gene_bulk_download", "Download Genes"),
                    )
                )
            ),
            # --- Tab: UMAP Generation ---
            ui.nav_panel("UMAP generation",
                ui.page_fluid(
                    ui.br(),
                    ui.h3("Generate UMAP by Gene"),
                    ui.p("Enter one or two genes to visuzalize their expression on a UMAP plot.", style="margin-bottom: 0.25rem;"),
                    ui.p("Note: plot genertion takes 1 to 2 minutes", style="margin-top: 0rem;"),
                    ui.br(),
                    ui.input_text("gene1", "Enter Gene 1", placeholder="EX: rab-3"),
                    ui.input_text("gene2", "Enter Gene 2 (optional for blend)"),
                    ui.input_switch(
                        "dataset_switch",
                        "Switch on for Male dataset",
                        value=False
                        ),
                    ui.input_action_button("run_umap","Generate UMAP")
                    )
            ),

            # --- Tab: heatmap generation ---
            ui.nav_panel("Heatmap generation",
                ui.page_fluid(
                    ui.br(),
                    ui.h3("Heatmap generation"),
                    ui.input_text_area("gene_names_heatmap", "Select all genes to include", placeholder="List each gene on its own line. You can also use an asterisk to get all genes with a prefix (ex. flp-*)", spellcheck="false", height="100px"),
                    ui.input_action_button("generate_heatmap", "Generate Heatmap"),
                )
            ),
            # -- Tab: Violin Plot Generation ---
            ui.nav_panel("Violin plot generation",
                    ui.page_fluid(
                        ui.br(),
                        ui.h3("Violin plot generation"),
                        ui.p("Plot generation takes 1 to 2 minutes", style="margin-bottom: 0.25rem;"),
                        
                        ui.input_text_area(
                            "gene_names",
                            "Enter gene(s) (one per line)" ,
                            spellcheck="false",
                            height="100px",
                            placeholder="e.g.\nrab-3\ntwk-18\nunc-17"
                            ),

                        ui.input_text_area(
                            "neuron_names",
                            "Enter Neuron(s) (one per line)",
                            spellcheck="false",
                            placeholder="e.g.\nAIY\nRIG\nASH", 
                            height="100px"
                            ),
                        
                        ui.input_switch(
                            "differential_switch",
                            "Switch on for differential comparison",
                            value=False
                            ),
                        ui.input_switch(
                            "violin_dataset_switch",
                            "Switch on for Male dataset",
                            value=False
                            ),

                        ui.p("Note: When showing differential expression it is N2 vs daf-2 or Male vs Hermaphrodite", style="margin-top: 0rem;"),

                        ui.input_action_button("run_violin", "Generate Violin Plot"),

                        ui.br(),
                        ui.output_ui("violin_output_ui")
                        )

        ),
                        

            # ---More Info Tab ---
            ui.nav_panel("More info",
                ui.page_fluid(
                    ui.br(),
                    ui.h3("Seurat Object"),
                    ui.p("The single nucleus data is an .rds file containing Seurat objects for all cells. This website contains many optimized queries for interfacing with the data, but if your analysis is more fine-grained, you may download the 7GB object below."),
                    ui.download_button("download_raw_data", "Download data"),

                    ui.div(
                        ui.h3("snSeq Papers by the Murphy Lab"),
                        ui.p("There are various circumstances where single nucleus RNA sequencing is the appropriate way to ask a biological question. "
                            "Here are some papers by our lab where we've employed this method."
                        ),
                        ui.layout_columns(
                            # Paper 1
                            ui.card(
                                ui.tags.a(
                                    "Adult single-nucleus neuronal transcriptomes of insulin signaling mutants reveal regulators of behavior and learning",
                                    href="https://doi.org/10.1016/j.xgen.2024.100720",
                                    target="_blank",
                                    style="font-weight: bold; font-size: 1.1rem;"
                                    ),
                                ui.tags.img(src="papers/StAnge_graphical_abstract.png", style="width: 100%; margin-top: 1rem; border-radius: 0.5rem;"),
                                ui.p("Gene expression in individual neurons can change during development to adulthood and can have large effects on behavior. Additionally, the insulin/insulin-like signaling (IIS) pathway regulates many of the adult functions of Caenorhabditis elegans, including learning and memory, via transcriptional changes. We used the deep resolution of single-nucleus RNA sequencing to define the adult transcriptome of each neuron in wild-type and daf-2 mutants, revealing expression differences between L4 larval and adult neurons in chemoreceptors, synaptic genes, and learning/memory genes. We used these data to identify adult new AWC-specific regulators of chemosensory function that emerge upon adulthood. daf-2 gene expression changes correlate with improved cognitive functions, particularly in the AWC sensory neuron that controls learning and associative memory; behavioral assays of AWC-specific daf-2 genes revealed their roles in cognitive function. Combining technology and functional validation, we identified conserved genes that function in specific adult neurons to control behavior, including learning and memory.")
                            ),

                            # Paper 2
                            ui.card(
                                ui.tags.a(
                                    "Body-to-brain insulin and Notch signaling regulates memory through neuronal CREB activity",
                                    href="https://doi.org/10.1038/s43587-025-00873-7",
                                    target="_blank",
                                    style="font-weight: bold; font-size: 1.1rem;"
                                ),
                                ui.tags.img(src="papers/Zhou_graphical_abstract.png", style="width: 100%; margin-top: 1rem; border-radius: 0.5rem;"),
                                ui.p("While memory regulation is predominantly understood as autonomous to neurons, factors outside the brain can also affect neuronal function. In Caenorhabditis elegans, the insulin/IGF-1-like signaling (IIS) pathway regulates longevity, metabolism and memory: long-lived daf-2 insulin/IGF-1 receptor mutants more than double memory duration after a single training session, and it was assumed that memory regulation was strictly neuronal. However, here we show that degradation of DAF-2 in the hypodermis also greatly extends memory, via expression of the diffusible Notch ligand, OSM-11, which in turn activates Notch signaling in neurons. Single-nucleus RNA sequencing of neurons revealed increased expression of CREB and other memory genes. Furthermore, in aged animals, activation of the hypodermal IIS–Notch pathway as well as OSM-11 overexpression rescue both memory and learning via CREB activity. Thus, insulin signaling in the liver-like hypodermis non-autonomously regulates neuronal function, providing a systemic connection between metabolism and memory through IIS–Notch–CREB signaling from the body to the brain.")
                            ),

                            # Paper 3
                            ui.card(
                                ui.tags.a(
                                    "Single-Nucleus Neuronal Transcriptional Profiling of Male C. elegans Uncovers Regulators of Sex-Specific and Sex-Shared Behaviors",
                                    href="https://doi.org/10.1016/j.celrep.2025.116016",
                                    target="_blank",
                                    style="font-weight: bold; font-size: 1.1rem;"
                                ),
                                ui.tags.img(src="papers/Morillo_graphical_abstract.png", style="width: 100%; margin-top: 1rem; border-radius: 0.5rem;"),
                                ui.p("Sexual differentiation of the nervous system causes differences in neuroanatomy, synaptic connectivity, and physiology. These sexually-dimorphic phenotypes ultimately translate into profound behavioral differences. C. elegans’ two sexes, XO males and XX hermaphrodites, demonstrate differences in neurobiology and behavior. However, the neuron class and sex-specific transcriptomic differences, particularly at the single-neuron level, that cause such phenotypic divergence remains understudied. Here, using single-nucleus RNA sequencing, we assessed and compared adult male and hermaphrodite C. elegans neuronal transcriptomes, identifying sex-specific neurons, including previously-unannotated male neurons. Sex-shared neurons displayed large expression differences, with some neuron classes clustering as distinct neurons between the sexes. Males express ∼100 male-specific GPCRs, largely limited to a subset of neurons. We identified the most highly-divergent neurons between the sexes, and functionally characterized a sex-shared target, vhp-1, in male-specific pheromone chemotaxis. Our data provide a resource for discovering nervous-system-wide sex transcriptomic differences and the molecular basis of sex-specific behaviors.")
                            ),

                            col_widths=[4, 4, 4],
                            style="margin-top: 2rem;"
                            ),
                        style="margin-top: 4rem; padding: 2rem; background-color: var(--bs-body-bg); border-radius: 1rem;"
                    )

                )
            ),
        ),

        #Output section (shown below the tabs)
        ui.div(
            ui.card(    
                ui.card_header("Output"),
                ui.row(
                    id="tables"
                ),
            ),
            style="margin-top: 2rem;"
        )
    ),
    title="Adult Worm Single Nucleus RNA Sequencing Atlas"
)

#TMP_IMAGE_PATH = "/tmp/user_umap.png"
#TMP_PNG_PATH = "/tmp/user_umap.png"


rendered_plots = {}

def server(input, output, session):
    mod_counter = reactive.value(0)

    @reactive.Effect
    @reactive.event(input.genotype)
    def update_cell_types():
        print("Updating cell types for genotype")
        new_choices = sorted([ct.upper() for ct in db.get_cell_types(input.genotype())])
        print(new_choices)

        ui.update_selectize(
            "cell_type",
            choices=new_choices,
        )

    @reactive.effect
    @reactive.event(input.add_query_gene_expr)
    def add_query_gene_expr():
        print("Running: add_query_gene_expr")
        selected_cell_type = input.cell_type()
        selected_threshold = input.threshold()
        results = db.get_gene_by_cell_type(input.genotype(), selected_threshold, selected_cell_type)
        counter = mod_counter.get() + 1
        mod_counter.set(counter)
        id = "query_" + str(counter)
        ui.insert_ui(
            selector="#tables",
            where="afterBegin",
            ui=ui.column(
                3,
                ui.tags.div(query_output_ui(id),
                    style="margin-bottom: 16px;"
                    ),
                id=id),
        )

        query_output_server(id, id, results, f"Cell Type={selected_cell_type}, {input.genotype()} threshold {selected_threshold}")
    
    @reactive.effect
    @reactive.event(input.add_query_cell_type)
    def add_query_cell_type():
        selected_cell_type = input.gene_name()
        selected_threshold = input.gene_threshold()
        results = db.get_cell_type_by_gene(input.genotype(), selected_threshold, selected_cell_type)
        counter = mod_counter.get() + 1
        mod_counter.set(counter)
        id = "query_" + str(counter)
        ui.insert_ui(
            selector="#tables",
            where="afterBegin",
            ui=ui.column(
                3,
                ui.tags.div(query_output_ui(id),
                    style="margin-bottom: 16px;"
                    ),
                id=id),
        )

        query_output_server(id, id, results, name=f"Gene={selected_cell_type}, {input.genotype()} threshold {selected_threshold}")
    
    @render.download(filename="bulk.csv")
    def gene_bulk_download():
        selected_genes = input.gene_bulk_names().splitlines()
        print(selected_genes)
        selected_threshold = input.gene_bulk_threshold()

        results = db.get_genes_bulk(input.genotype(), selected_threshold, selected_genes)
        print(results)
    
        if results is None:
            print("No results returned")
        
        yield db.df_to_csv(results)
    

    
   # @output
   # @render.ui
   # def user_umap_output():
    #    path = current_umap_path()
     #   if not path:
      #      return ui.div("No UMAP Plot yet")
       # print(f"Rendering umap plot for path: {path}")
       # return ui.tags.img(src=path)

    @reactive.effect
    @reactive.event(input.run_umap)
    def display_umap():
        gene1 = input.gene1()
        gene2 = input.gene2()
        data_toggle = input.dataset_switch()

        #Check if gene 1 has been inputted by user
        if not gene1 or not gene1.strip():
            print("Gene 1 is required")
            return

        gene1 = gene1.strip()
        gene2 = gene2.strip()
        
        #Generate a unique PNG filename in the www folder
        # Generate a unique ID and path
        counter = mod_counter.get() + 1
        mod_counter.set(counter)
        query_id = f"query_{counter}"
        unique_filename = f"umap_{uuid.uuid4().hex}.png"
        output_path = os.path.join(CESTAAN_ROOT, "www", unique_filename)
        

        #Call R script with 1 or 2 genes and output path
        r_script = os.path.join(CESTAAN_ROOT, "generate_umap.R")
        cmd = ["Rscript", r_script, gene1, output_path]

        if gene2:
            cmd.extend(["--gene2", gene2])
        
        if data_toggle:
            cmd.append("--dataset")


        print(f"Running command: {' '.join(cmd)}")

        try:
            subprocess.run(cmd, check=True)
            print("UMAP image generated successfully.")
            print(f"[DEBUG] Just generated image: {output_path}")
            print(f"[DEBUG] Setting current_umap_path to: /{unique_filename}")
            
        except FileNotFoundError:
            print("Error: R script not found or not executable.")
            return        
        except subprocess.CalledProcessError as e:
            print(f"Rscript failed:{e}")
            return
        
       # relative_web_path = "/" + unique_filename
       # current_umap_path.set(relative_web_path)

        #counter = mod_counter.get() + 1
        #mod_counter.set(counter)
        #query_id = f"query_{counter}"
        
        ui.insert_ui(
            selector="#tables",
            where="afterBegin",
            ui=ui.div(
                query_umap_ui(query_id),
                id=query_id 
                )
            )
        
        img_path = unique_filename
        #print(f"[INFO 1A] Inserting query {query_id} - output function is: {user_umap_output}")
        print(f"[INFO 1B] query+umap_ui({query_id}) returned {query_umap_ui(query_id)}")
       # print(f"[INFO 1C] name of output path to call image: {relative_web_path}")
        #query_umap_server(query_id, query_id, f"/www/{current_umap_path()}", name="UMAP Plot:") 
        print(f"[INFO 1D] unique filename: {img_path}")
        query_umap_server(
                query_id, 
                remove_id = query_id, 
                img_path = img_path, 
                name = f"UMAP Plot: {gene1}"+ (f" + {gene2}" if gene2 else ""), 
                is_two_genes=bool(gene2)
            )

    print("Look here for potential error")
    print(f"UMAP UI id: {id}, type: {type(id)}")
    
    # Commented out used for testing reactive displays
    #@reactive.effect
    #@reactive.event(input.test_display)
    #def test_display_umap():
    #    counter = mod_counter.get() + 1
    #    mod_counter.set(counter)
    #    query_id = f"test_query_{counter}"
    #    print(query_id)

    #    ui.insert_ui(
    #            selector= "#tables",
    #            where = "beforeEnd",
    #            ui= query_umap_ui(query_id)
    #            )
    #    query_umap_server(query_id, remove_id = query_id, img_path = TEST_IMAGE_PATH, name = "Static Test UMAP")
   
    @reactive.effect
    @reactive.event(input.generate_heatmap)
    def generate_heatmap():
        print("generating heatmap")
        selected_genes = input.gene_names_heatmap().splitlines()

        if not selected_genes:
            print("No genes entered.")
            return
        
        avg, percent = db.get_avg_percent_for_genes(input.genotype(), selected_genes)

        plot = heatmap.generate_interactive_heatmap(avg, percent)
        #fig = heatmap.generate_interactive_heatmap(avg, percent)


        counter = mod_counter.get() + 1
        mod_counter.set(counter)
        id = "query_" + str(counter)
        ui.insert_ui(
            selector="#tables",
            where="afterBegin",
            ui=ui.column(
                12,
                ui.tags.div(query_heatmap_ui(id),
                    style="margin-bottom: 16px; padding: 8px; boarder: 2px solid #ccc;"
                    ),
                id=id),
        )

        query_heatmap_server(id, id, plot, name=f"Plot from {input.genotype()}")
    
    @reactive.effect
    @reactive.event(input.run_violin)
    def display_violin():
        genes_raw = input.gene_names()
        neurons_raw = input.neuron_names()
        diff_toggle = input.differential_switch()
        dataset_toggle = input.violin_dataset_switch()

        if not genes_raw or not neurons_raw:
            print("Both genes and neurons must be provided")
            return

    # Clean up inputs
        genes = ",".join([g.strip() for g in genes_raw.splitlines() if g.strip()])
        neurons = ",".join([n.strip() for n in neurons_raw.splitlines() if n.strip()])

        if not genes or not neurons:
            print("Invalid input: genes or neurons list is empty after cleaning")
            return
                                 
        # Generate unique output path
        counter = mod_counter.get() + 1
        mod_counter.set(counter)
        
        query_id = f"query_violin_{counter}"
        filename = f"violin_{uuid.uuid4().hex}.png"
        
        output_rel_path = f"www/{filename}"
        output_abs_path = os.path.join(CESTAAN_ROOT, output_rel_path)

        # Build Rscript command
        #r_script = "/var/www/ctvm1/murphy-lab-project/generate_violin.R"
        #cmd = ["Rscript", r_script, genes, neurons, output_abs_path]
        
        # Construct R command
        r_command = [
                "Rscript",
                "generate_violin.R",
                "--genes", genes,
                "--neurons", neurons,
                "--output", output_abs_path,
            ]

        #print(f"Running command: {' '.join(cmd)}")
                                 

        if diff_toggle:
            r_command.append("--differential")
        if dataset_toggle:
            r_command.append("--dataset")
        
        print(f"[INFO] Running R command: {' '.join(r_command)}")

            
        try:
            result = subprocess.run(r_command,capture_output=True, text=True, check=True)
            print("Violin plot generated successfully.")
            print(f"[DEBUG] Violin plot image at: {output_abs_path}")
        except FileNotFoundError:
            print("Error: R script not found or not executable.")
            return
        except subprocess.CalledProcessError as e:
            print(f"Rscript failed: {e}")
            return
        
        ui.insert_ui(
                selector="#tables",
                where="afterBegin",
                ui=ui.div(
                    query_umap_ui(query_id),
                    id=query_id
                    )
            )

        query_umap_server(
                id=query_id,
                remove_id=query_id,
                img_path=output_rel_path,
                name=f"Violin Plot: {genes} in {neurons}",
                is_two_genes=(len(genes.split(",")) > 1)
            )
  


    @render.download()
    def download_raw_data():
        return "source_data/WT_daf2.rds"
    
    @render.image
    def umap_image():
        img = {"src": here / "umap.png", "height": "450px"}  
        return img

#from fastapi.staticfiles import StaticFiles
#app.mount("/files", StaticFiles(directory="/tmp/violin_output"), name="files")

db.init()
app = App(app_ui, server, static_assets= static_dir)

