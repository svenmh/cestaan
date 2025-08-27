import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
import numpy as np

def generate_interactive_heatmap(
    avg_df: pd.DataFrame, 
    prop_df: pd.DataFrame
):
    """
    Generates an interactive heatmap using Plotly from the given DataFrames.
    
    Parameters:
    - avg_df: pd.DataFrame : DataFrame containing the average expression levels (e.g., TPM).
    - prop_df: pd.DataFrame : DataFrame containing the proportions (percentages).
    
    Assumes:
    - The first column is 'gene'.
    - The remaining columns in avg_df and prop_df correspond to cell types.
    
    Returns:
    - fig: plotly.graph_objects.Figure : The Plotly figure object for the heatmap.
    """
    
    # Melt the DataFrames to long format
    print("avg_df shape:", avg_df.shape)
    print("avg_df columns:", avg_df.columns)
    print("avg_df preview:")
    print(avg_df.head())
    avg_df_melted = avg_df.melt(id_vars=[avg_df.columns[0]], var_name='CellType', value_name='Expression')
    print("prop_df shape:", prop_df.shape)
    print("prop_df columns:", prop_df.columns)
    print("prop_df preview:")
    print(prop_df.head())
    prop_df_melted = prop_df.melt(id_vars=[prop_df.columns[0]], var_name='CellType', value_name='Proportion')
    
    # Merge the two DataFrames on the gene and cell type columns
    merged_df = pd.merge(avg_df_melted, prop_df_melted, on=[avg_df.columns[0], 'CellType'])
    print(type(merged_df))
    print(list(merged_df['CellType']))
    
    #Get the max proportion to scale dotsizes properly (this will help with legend)
    max_prop = merged_df["Proportion"].max()
    #These are based on how plotly scales dotsizes
    desired_max_size = 10
    size_ref = max_prop / desired_max_size**2

    # Create the scatter plot with adjusted parameters to avoid overlap
    fig = px.scatter(
        merged_df, 
        x='CellType', 
        y=avg_df.columns[0],  # The 'gene' column
        size='Proportion', 
        color='Expression',
        hover_name=avg_df.columns[0],  # What you see on hover (gene)
        hover_data={'Expression': True, 'Proportion': True},  # Additional hover data
        color_continuous_scale='Plasma_r',  # You can choose another color scale
        size_max=desired_max_size,  # Reduce max size to prevent overlap
        height = max(400, 0.2 * len(avg_df_melted.index) + 5)
    )
    fig.update_traces(marker=dict(sizeref=size_ref,sizemode='area'))

    # --- Add dot size legend (based on actual Proportion values) ---
    
# Define example points for legend using percentiles or fixed values
   # legend_props = np.percentile(merged_df['Proportion'], [10, 50, 90])
    legend_props = [
            merged_df['Proportion'].min(),
            merged_df['Proportion'].mean(),
            merged_df['Proportion'].max()
            ]
    legend_labels = [f"{int(p)}%" for p in legend_props]

# px.scatter uses size_max=10 and automatically rescales proportions to [0, size_max]
# so we mimic that here

    def scale_size(prop):
        if legend_props[0] == legend_props[2]:
            return desired_max_size / 2  # fallback if all values are equal
        return (prop/legend_props[2])*desired_max_size

    for prop, label in zip(legend_props, legend_labels):
        fig.add_trace(go.Scatter(x=[None], y=[None],
            mode='markers',
            marker=dict(
                size=scale_size(prop),
                color='gray',
                opacity=0.6) ,
                showlegend=True,
                name=f"{label} expressing"
        ))
                                                                                            

    # Adjust the layout
    fig.update_layout(
        xaxis=dict(
            tickfont=dict(size=12),  # Reduce font size to half the default (which is usually 16)
            title=dict(text="Cell Type",font=dict(size=14)),
            automargin=True,  # Ensure labels fit without rotation
            tickangle=90,  # Set tick angle to 90 for rotation
        ),
        yaxis=dict(
            tickfont=dict(size=12),  # Reduce font size to half the default (which is usually 16)
            title=dict(text="Gene", font=dict(size=14)),
        ),
        legend=dict(
            title="Proportion Expressing",
           # orientation = "h", #horizontal legend
            x=1, #Left and right
            y=-0.5,
            traceorder="normal",
            itemsizing='constant'
        ),
        coloraxis_colorbar=dict(
            title="Expression",
            x=1 #Left and right
        )
    )

    # Return the figure object
    return fig
