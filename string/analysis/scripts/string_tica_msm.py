import sys

import matplotlib.pyplot as plt
import numpy as np
from stringmethod.config import *
from stringmethod.postprocessing import *

plt.rcParams["axes.facecolor"] = "#f9f9fb"
plt.rcParams["grid.color"] = "white"
plt.rcParams["grid.linestyle"] = "-"
plt.rcParams["grid.linewidth"] = 2
plt.rcParams["axes.grid"] = True
plt.rcParams["lines.solid_capstyle"] = "round"


def _colorbar(mappable, cmap, norm, label0, size=10):
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    ax = mappable.axes
    # fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cbar.set_label(label0, size=size)
    return cbar


def load_swarm_data(extract, first_iteration=1, last_iteration=None):
    if last_iteration is None:
        last_iteration = sys.maxsize
    if extract:
        config = load_config("config_postprocessing.json")

        ce = CvValueExtractor.from_config(
            config=config,
            # Exclude the first iterations to let the system equilibrate.
            first_iteration=first_iteration,
            # Usefull to make blocks of the simulation
            last_iteration=last_iteration,
        )
        ce.run()
        ce.persist()
    return np.load("postprocessing/cv_coordinates.npy")


def cvs_to_tica(cv_coordinates, drop):
    from deeptime.decomposition import TICA

    data = []
    cvs = list([i for i in range(cv_coordinates.shape[2]) if i not in drop])
    for i in range(cv_coordinates.shape[0]):
        data.append(cv_coordinates[i, :, cvs].T)

    tica = TICA(lagtime=1)
    data = tica.fit(data, lagtime=1, progress=True).fetch_model().transform(data)
    return data


def k_means_cluster(data, k, stride=1, max_iter=500, n_proc=1, seed=None):
    from deeptime.clustering import KMeans

    estimator = KMeans(
        n_clusters=k,  # place 100 cluster centers
        init_strategy="uniform",  # uniform initialization strategy
        max_iter=max_iter,  # don't actually perform the optimization, just place centers
        fixed_seed=seed,
        n_jobs=n_proc,
    )
    clusters = (
        estimator.fit(data[::stride, :, :].reshape(-1, data.shape[2]))
        .fetch_model()
        .transform(data.reshape(-1, data.shape[2]))
    )
    clusters = clusters.reshape(-1, data.shape[1])
    return clusters


def get_vamp_vs_k(n_clustercenters, data):
    import logging

    import deeptime.markov as markov
    from deeptime.decomposition import vamp_score_cv
    from deeptime.util import confidence_interval
    from tqdm.autonotebook import tqdm

    loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
    for logger in loggers:
        logger.setLevel(logging.ERROR)

    n_iter = 5
    scores = np.zeros((len(n_clustercenters), n_iter))
    for n, k in tqdm(
        enumerate(n_clustercenters), total=len(n_clustercenters), desc="Loop over k:"
    ):
        for m in tqdm(range(n_iter), desc="Loop over iterations:", leave=False):
            _cl = k_means_cluster(data, k, stride=10, max_iter=50, n_proc=8)

            estimator = markov.msm.MaximumLikelihoodMSM(
                reversible=True,
                stationary_distribution_constraint=None,
                lagtime=1,
            )

            counts = (
                markov.TransitionCountEstimator(lagtime=1, count_mode="sample")
                .fit(_cl)
                .fetch_model()
            )
            _msm = estimator.fit(counts)
            # return _msm, _cl
            # exit

            scores[n, m] = vamp_score_cv(
                _msm, trajs=[c for c in _cl], n=1, lagtime=1, dim=min(10, k)
            )[0]

        # Plotting
    fig, ax = plt.subplots(1, 1)
    lower, upper = confidence_interval(scores.T.tolist(), conf=0.9)
    ax.fill_between(n_clustercenters, lower, upper, alpha=0.3)
    ax.plot(n_clustercenters, np.mean(scores, axis=1), "-o")
    ax.semilogx()
    ax.set_xlabel("number of cluster centers")
    ax.set_ylabel("VAMP-2 score")
    fig.tight_layout()

    return fig, ax


def get_bayesian_msm(clusters, n_samples=100):
    import deeptime.markov as markov

    estimator = markov.msm.BayesianMSM(
        n_samples=100, reversible=True, stationary_distribution_constraint=None
    )

    counts = (
        markov.TransitionCountEstimator(lagtime=1, count_mode="effective")
        .fit(clusters)
        .fetch_model()
    )

    msm = estimator.fit(counts).fetch_model()
    print("here")
    weights = msm.gather_stats(
        "compute_trajectory_weights", dtrajs=clusters.reshape(-1, 1)
    )
    print("here")

    return msm, weights


def get_msm(clusters):
    import logging

    import deeptime.markov as markov

    loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
    for logger in loggers:
        logger.setLevel(logging.ERROR)

    estimator = markov.msm.MaximumLikelihoodMSM(
        reversible=True,
        stationary_distribution_constraint=None,
        lagtime=1,
    )
    counts = (
        markov.TransitionCountEstimator(lagtime=1, count_mode="effective")
        .fit(clusters)
        .fetch_model()
    )

    msm = estimator.fit(counts).fetch_model()
    weights = np.array(
        msm.compute_trajectory_weights(np.concatenate(clusters.reshape(-1, 1)))
    ).squeeze()

    return msm, weights


def get_kde(samples, weights=None, bandwidth=None, nbins=55, extent=None):
    """
    aklsjfa;sldjfas
    """

    from scipy.stats import gaussian_kde

    samples = samples.reshape(-1, samples.shape[1])
    if extent is None:
        xmin = samples[:, 0].min()
        xmax = samples[:, 0].max()
        ymin = samples[:, 1].min()
        ymax = samples[:, 1].max()
    else:
        xmin, xmax, ymin, ymax = extent
    nbins = nbins * 1j
    X, Y = np.mgrid[xmin:xmax:nbins, ymin:ymax:nbins]
    positions = np.vstack([X.ravel(), Y.ravel()])
    kernel = gaussian_kde(
        samples.T,
        weights=weights,
        bw_method=bandwidth,
    )
    Z = np.reshape(kernel(positions), X.shape)
    Z = Z.T
    extent = [xmin, xmax, ymin, ymax]

    return Z, extent


def get_error(
    cv_proj,
    clusters,
    extent,
    n_boot=200,
    bandwidth=0.05,
    nbin=55,
    blocks=[2, 4, 8, 16, 32, 64],
):
    from tqdm.autonotebook import tqdm

    assert (
        type(bandwidth) == float
    ), "You have to choose a bandwidth, otherwise the bandwith will change between bootstrap iterations"
    ndat = cv_proj.shape[0]
    errors = []
    for b in tqdm(blocks, desc="Loop over blocks"):
        histograms = []
        block_length = ndat // b
        for _ in tqdm(range(n_boot), leave=False, desc="Loop over bootstraps"):
            random = np.random.choice(b, b)
            mask = []
            for r in random:
                mask += list(range(r * block_length, (r + 1) * block_length))
            _, w = get_msm(clusters[mask])
            h, extent = get_kde(
                cv_proj[mask, :, :],
                w,
                bandwidth,
                extent=extent,
                nbins=nbin,
            )
            histograms.append(h)
        histograms = np.array(histograms)
        x_mean, hdi = get_hdi(histograms, 0, 0.05)
        # Dividing by x_mean propagates the uncertainty from histogram uncertainty to
        # free energy uncertainty.
        errors.append((hdi[1, :, :] - hdi[0, :, :]) / 2 / x_mean)
    errors = np.array(errors)
    errors[errors == np.inf] = np.nan

    return errors


def get_hdi(x, axis, alpha=0.06):
    x_mean = np.nanmean(x, axis=axis)
    percentiles = 100 * np.array([alpha / 2.0, 1.0 - alpha / 2.0])
    hdi = np.nanpercentile(x, percentiles, axis=axis)
    return x_mean, hdi


def plot_2D_heatmap(
    G,
    extent,
    cmap=None,
    f_min=0,
    f_max=20,
    cbar_label="Free Energy(kT)",
    ax=None,
    fig=None,
    xlabel="",
    ylabel="",
    xlim=None,
    ylim=None,
):
    import matplotlib as mpl

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 7), sharex=True, sharey=True)
    if cmap is None:
        cmap = plt.cm.RdYlBu_r
    n_colors = 50
    colors = cmap(np.linspace(0, 1, n_colors))  # yellow to blue
    norm = mpl.colors.Normalize(vmin=f_min, vmax=f_max)
    ax.contourf(
        G,
        cmap=cmap,
        extent=extent,
        vmax=f_max,
        vmin=f_min,
        levels=n_colors,
    )
    _ = _colorbar(ax, cmap, norm, cbar_label, 15)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    fig.tight_layout()
    return fig, ax
