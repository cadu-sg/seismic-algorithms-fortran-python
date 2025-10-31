import colorcet as cc
import numpy as np
import numpy.typing as npt
from bokeh.layouts import row
from bokeh.models import ColumnDataSource, GlyphRenderer, PointDrawTool
from bokeh.plotting import figure, show, curdoc
from seismicio import readsu

import seismic_algorithms

indatafile = "/storage1/Seismic/dados_teste/marmousi_4ms_CDP.su"

sufile = readsu(indatafile, gather_keyword="cdp")

# PREPARE DATA
# ------------

# Number of time samples per data trace
num_time_samples = sufile.num_samples

# Number of traces
num_traces = sufile.gather[100].num_traces

cdp_gather_data = sufile.gather[100].data
cdp_gather_offsets = sufile.gather[100].headers["offset"]

first_time_sample = 0.0
interval_time_samples = (
    sufile.headers["dt"][0] / 1e6
)  # convert microseconds (µs) to seconds (s)
last_time_sample = first_time_sample + (num_time_samples - 1) * interval_time_samples
width_time_samples = np.abs(last_time_sample - first_time_sample)

width_offsets = np.abs(cdp_gather_offsets[0] - cdp_gather_offsets[-1])

# CREATE CDP GATHER PLOT
# ----------------------

# Create figure object
cdp_gather_plot = figure(
    active_drag=None,
    x_axis_label="Offset (m)",
    x_axis_location="above",
    y_axis_label="Time (s)",
)

# Adjust plot ranges
cdp_gather_plot.x_range.range_padding = 0.0
cdp_gather_plot.y_range.range_padding = 0.0
cdp_gather_plot.y_range.flipped = True

cdp_gather_renderer = cdp_gather_plot.image(
    image=[cdp_gather_data],
    x=cdp_gather_offsets[0],
    y=first_time_sample,
    dw=width_offsets,
    dh=width_time_samples,
    anchor="bottom_left",
    origin="bottom_left",
    palette="Greys256",
)

# %%

# %%
# COMPUTE SEMBLANCE DATA
# ----------------------


def semblance(
    cdp_gather_data: npt.NDArray,
    cdp_gather_offsets: npt.NDArray,
    velocities: npt.NDArray,
    first_time_sample: float,
    interval_time_samples: float,
) -> npt.NDArray:
    num_time_samples = cdp_gather_data.shape[0]
    num_traces = cdp_gather_data.shape[1]

    return seismic_algorithms.seismic_algorithms.semblance(
        sucmpdata=cdp_gather_data,
        offsets=cdp_gather_offsets,
        velcoer=velocities,
        t0_data=first_time_sample,
        dt=interval_time_samples,
        nt=num_time_samples,
        ntraces=num_traces,
        nvelcoer=len(velocities),
    )


vel_min = 1000.0
vel_max = 5000.0
vel_step = 25.0
velocities = np.arange(vel_min, vel_max + 0.1, vel_step, dtype=float)
width_velocities = np.abs(vel_max - vel_min)

coherence_matrix = semblance(
    cdp_gather_data=cdp_gather_data,
    cdp_gather_offsets=cdp_gather_offsets,
    velocities=velocities,
    first_time_sample=first_time_sample,
    interval_time_samples=interval_time_samples,
)

print(f"min: {np.min(coherence_matrix)}, max: {np.max(coherence_matrix)}")

# CREATE SEMBLANCE PLOT
# ---------------------

# Create figure object
semblance_plot = figure(
    active_drag=None,
    x_axis_label="Velocities (m/s)",
    x_axis_location="above",
    y_axis_label="Time (s)",
)

# Adjust plot ranges
semblance_plot.x_range.range_padding = 0.0
semblance_plot.y_range.range_padding = 0.0
semblance_plot.y_range.flipped = True

semblance_renderer: GlyphRenderer = semblance_plot.image(
    image=[coherence_matrix],
    x=velocities[0],
    y=first_time_sample,
    dw=width_velocities,
    dh=width_time_samples,
    anchor="bottom_left",
    origin="bottom_left",
    palette=cc.rainbow4,
)

semblance_color_bar = semblance_renderer.construct_color_bar(padding=1)

semblance_plot.add_layout(semblance_color_bar, "right")

# POINT DRAWER
# ------------

semblance_picks_source = ColumnDataSource(data={"x": [], "y": []})

semblance_picks_scatter_renderer = semblance_plot.scatter(
    x="x",
    y="y",
    color="black",
    source=semblance_picks_source,
)

semblance_plot.line(
    x="x",
    y="y",
    color="black",
    source=semblance_picks_source,
    line_width=1.5,
)

semblance_plot.add_tools(PointDrawTool(renderers=[semblance_picks_scatter_renderer]))


def sort_xy_pairs_by_x(x, y):
    x_array = np.asarray(x)
    y_array = np.asarray(y)
    sorted_indices = np.argsort(x_array)
    return x_array[sorted_indices], y_array[sorted_indices]



is_updating_picks_data = False
def picks_data_on_change_handler(attr, old, new):
    global is_updating_picks_data
    if is_updating_picks_data:
        return
    is_updating_picks_data = True
    x_sorted, y_sorted = sort_xy_pairs_by_x(new["x"], new["y"])
    semblance_picks_source.data = {"x": x_sorted, "y": y_sorted}
    is_updating_picks_data = False


semblance_picks_source.on_change("data", picks_data_on_change_handler)

# SHOW PLOTS SIDE BY SIDE
# -----------------------

row_layout_plots = row(children=[cdp_gather_plot, semblance_plot])

# show(row_layout_plots)

curdoc().add_root(row_layout_plots)
