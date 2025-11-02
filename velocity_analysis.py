import colorcet as cc
import numpy as np
import numpy.typing as npt
from bokeh.layouts import row, column
from bokeh.models import ColumnDataSource, GlyphRenderer, PointDrawTool, Button
from bokeh.plotting import figure, curdoc
from seismicio import readsu

import seismic_algorithms


# VELOCITY ANALYSIS ALGORITHMS WRAPPERS
# -------------------------------------


def velocity_picks_to_trace(
    picks_times: npt.ArrayLike,
    picks_velocities: npt.ArrayLike,
    first_time_sample: float,
    interval_time_samples: float,
    cdp_gather_num_time_samples: int,
) -> npt.NDArray:
    return seismic_algorithms.seismic_algorithms.velpicks_to_trace(
        npicks=len(picks_times),
        tnmo=picks_times,
        vnmo=picks_velocities,
        t0_data=first_time_sample,
        dt=interval_time_samples,
        nt=cdp_gather_num_time_samples,
    )


def apply_nmo(
    cdp_gather_data: npt.NDArray,
    offsets: npt.NDArray,
    interpolated_velocities_trace: npt.NDArray,
    first_time_sample: float,
    interval_time_samples: float,
    smute: float,
) -> npt.NDArray:
    return seismic_algorithms.seismic_algorithms.apply_nmo(
        ntracescmp=cdp_gather_data.shape[1],  # Number of traces of the CMP gather
        nt=cdp_gather_data.shape[0],  # Number of time samples of the CMP gather
        t0_data=first_time_sample,
        dt=interval_time_samples,
        cmpdata=cdp_gather_data,
        offsets=offsets,
        vnmo_trace=interpolated_velocities_trace,
        smute=smute,
    )


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


# PREPARE DATA
# ------------

indatafile = "/storage1/Seismic/dados_teste/marmousi_4ms_CDP.su"

sufile = readsu(indatafile, gather_keyword="cdp")

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

# COMPUTE SEMBLANCE DATA
# ----------------------


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


# CREATE NMO PLOT
# ---------------

# Create figure object
nmo_plot = figure(
    active_drag=None,
    x_axis_label="Offset (m)",
    x_axis_location="above",
    y_axis_label="Time (s)",
)

# Adjust plot ranges
nmo_plot.x_range.range_padding = 0.0
nmo_plot.y_range.range_padding = 0.0
nmo_plot.y_range.flipped = True

nmo_renderer = nmo_plot.image(
    image=[cdp_gather_data],
    x=cdp_gather_offsets[0],
    y=first_time_sample,
    dw=width_offsets,
    dh=width_time_samples,
    anchor="bottom_left",
    origin="bottom_left",
    palette="Greys256",
)

nmo_source : ColumnDataSource = nmo_renderer.data_source
print(type(nmo_source))
print(type(nmo_source.data))
print(type(nmo_source.data['image']))

# POINT DRAWER
# ------------

semblance_picks_source = ColumnDataSource(data={"x": [], "y": []})

semblance_picks_scatter_renderer = semblance_plot.scatter(
    x="x",
    y="y",
    color="black",
    source=semblance_picks_source,
)

semblace_picks_line_renderer = semblance_plot.line(
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


# is_updating_picks_data = False
def picks_data_on_change_handler(attr, old, new):
    # print(f"called {datetime.datetime.now()}")
    semblance_picks_source.remove_on_change("data", picks_data_on_change_handler)
    x_sorted, y_sorted = sort_xy_pairs_by_x(new["x"], new["y"])

    semblance_picks_source.data.update({"x": x_sorted, "y": y_sorted})

    semblance_picks_source.on_change("data", picks_data_on_change_handler)


semblance_picks_source.on_change("data", picks_data_on_change_handler)


def test_button_handler(new):
    print(semblance_picks_source.data)
    # print(interval_time_samples)
    interpolated_velocities_trace = velocity_picks_to_trace(
        picks_times=semblance_picks_source.data["y"],
        picks_velocities=semblance_picks_source.data["x"],
        first_time_sample=first_time_sample,
        interval_time_samples=interval_time_samples,
        cdp_gather_num_time_samples=num_time_samples,
    )

    print(type(interpolated_velocities_trace))

    nmo_corrected_cdp_gather_data = apply_nmo(
        cdp_gather_data=cdp_gather_data,
        offsets=cdp_gather_offsets,
        interpolated_velocities_trace=interpolated_velocities_trace,
        first_time_sample=first_time_sample,
        interval_time_samples=interval_time_samples,
        smute=1.5,
    )

    nmo_source.data = {"image": [nmo_corrected_cdp_gather_data]}


button = Button(label="APPLY NMO CORRECTION")
button.on_click(test_button_handler)


# ORGANIZE LAYOUT
# ---------------

row_layout_plots = row(children=[cdp_gather_plot, semblance_plot, nmo_plot])

outer_layout = column(row_layout_plots, button)


# TEST BUTTON


# show(row_layout_plots)

curdoc().add_root(outer_layout)
