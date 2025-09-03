import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

def plot_clustermap(infolder, config=None):

    OUTPUT_POCKETS_CSV = os.path.join(infolder, "pockets_clustered.csv")

    df_pockets_clustered = pd.read_csv(OUTPUT_POCKETS_CSV)
    
    # Exclude cluster -1 from the analysis (that is the "noise" cluster)
    df_pockets_clustered = df_pockets_clustered[df_pockets_clustered['cluster'] != -1]
    
    # --- Parse residues for each frame+cluster ---
    def parse_residues(res_str):
        if pd.isnull(res_str):
            return set()
        return set([r for r in str(res_str).split() if r])

    df_pockets_clustered['residue_set'] = df_pockets_clustered['residues'].apply(parse_residues)

    # Get all unique residues in this simulation
    all_residues = sorted(set.union(*df_pockets_clustered['residue_set']) if len(df_pockets_clustered) > 0 else set())

    # Sort by cluster, then by frame
    df_pockets_clustered.sort_values(['cluster', 'Frame'], inplace=True)

    unique_clusters = df_pockets_clustered['cluster'].unique()
    n_clusters = len(unique_clusters)

    # Calculate number of rows (frames) per cluster
    cluster_sizes = []
    cluster_residues = {}  # Store residues for each cluster
    for clust in unique_clusters:
        cluster_data = df_pockets_clustered[df_pockets_clustered['cluster'] == clust]
        cluster_sizes.append(len(cluster_data))
        # Get unique residues for this cluster
        cluster_residues[clust] = sorted(set.union(*cluster_data['residue_set']) if len(cluster_data) > 0 else set())
    
    cluster_sizes = np.array(cluster_sizes)
    # Use ratios for height allocation, but cap total height at 14
    min_height = 0.5  # Minimum height per cluster
    total_height = 14
    # Calculate raw ratios
    raw_ratios = cluster_sizes / cluster_sizes.sum()
    # Initial heights
    heights = raw_ratios * total_height
    # Enforce minimum height
    heights = np.maximum(heights, min_height)
    # Renormalize if sum exceeds total_height due to min_height enforcement
    if heights.sum() > total_height:
        heights = heights * (total_height / heights.sum())

    fig, axes = plt.subplots(
        nrows=n_clusters, ncols=1,
        figsize=(18, total_height),
        sharex=True, sharey=False,
        gridspec_kw={'hspace': 0, 'wspace': 0, 'height_ratios': heights}
    )
    if n_clusters == 1:
        axes = [axes]
    cmap_list = sns.color_palette("husl", n_clusters)

    # Store cluster boundaries and data for interaction
    cluster_boundaries = {}
    current_y = 0
    
    for idx, clust in enumerate(unique_clusters):
        ax = axes[idx]
        clust_df = df_pockets_clustered[df_pockets_clustered['cluster'] == clust]
        index = clust_df['Frame']
        onehot = pd.DataFrame(
            [
                [1 if res in row else 0 for res in all_residues]
                for row in clust_df['residue_set']
            ],
            index=index,
            columns=all_residues
        )
        
        # Store cluster boundary information
        cluster_boundaries[clust] = {
            'ax': ax,
            'y_start': current_y,
            'y_end': current_y + len(clust_df),
            'residues': cluster_residues[clust],
            'onehot': onehot
        }
        current_y += len(clust_df)
        
        sns.heatmap(
            onehot,
            cmap=sns.light_palette(cmap_list[idx], as_cmap=True),
            cbar=False,
            ax=ax,
            yticklabels=False,
            xticklabels=False,  # We'll handle x-labels interactively
        )
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title("")  # No axis titles

        # Add cluster label to the left of each subplot
        ax.text(
            -0.02, 0.5, f"Cluster {clust}",
            va='center', ha='right',
            fontsize=14, fontweight='bold',
            transform=ax.transAxes,
            rotation=0
        )

        ax.margins(0)

    # Interactive functionality
    class InteractiveClusterMap:
        def __init__(self, fig, axes, cluster_boundaries, all_residues):
            self.fig = fig
            self.axes = axes
            self.cluster_boundaries = cluster_boundaries
            self.all_residues = all_residues
            self.tooltip = None
            
            # Connect events
            self.fig.canvas.mpl_connect('motion_notify_event', self.on_hover)
            self.fig.canvas.mpl_connect('axes_leave_event', self.on_leave)
            
        def on_hover(self, event):
            if event.inaxes is None:
                self.clear_tooltip()
                return
                
            # Find which cluster we're hovering over
            for cluster_id, info in self.cluster_boundaries.items():
                if event.inaxes == info['ax']:
                    # Get the heatmap data for this cluster
                    onehot = info['onehot']
                    
                    # Convert mouse coordinates to data coordinates
                    # event.xdata and event.ydata are in data coordinates
                    if event.xdata is not None and event.ydata is not None:
                        x_idx = int(event.xdata)
                        y_idx = int(event.ydata)
                        
                        # Check bounds
                        if 0 <= x_idx < len(self.all_residues) and 0 <= y_idx < len(onehot):
                            residue = self.all_residues[x_idx]
                            
                            # Check if this residue is present in this position (value = 1)
                            if onehot.iloc[y_idx, x_idx] == 1:
                                frame_name = onehot.index[y_idx]
                                self.show_tooltip(event, residue, cluster_id, frame_name)
                            else:
                                self.clear_tooltip()
                        else:
                            self.clear_tooltip()
                    break
            else:
                self.clear_tooltip()
                    
        def on_leave(self, event):
            self.clear_tooltip()
            
        def show_tooltip(self, event, residue, cluster_id, frame_name):
            # Clear existing tooltip
            self.clear_tooltip()
            
            # Create tooltip text
            tooltip_text = f"Residue: {residue}\nCluster: {cluster_id}\nFrame: {frame_name}"
            
            # Create tooltip annotation
            self.tooltip = event.inaxes.annotate(
                tooltip_text,
                xy=(event.xdata, event.ydata),
                xytext=(20, 20),  # Offset from mouse position
                textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.8),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'),
                fontsize=10,
                zorder=1000
            )
            
            self.fig.canvas.draw_idle()
            
        def clear_tooltip(self):
            if self.tooltip:
                self.tooltip.remove()
                self.tooltip = None
                self.fig.canvas.draw_idle()

    # Create interactive handler
    interactive_handler = InteractiveClusterMap(fig, axes, cluster_boundaries, all_residues)

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.suptitle(f"Residue Composition per Cluster (Hover over heatmap points to see residue details)", y=1.01)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.show()