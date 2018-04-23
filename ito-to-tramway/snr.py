# -*- coding: utf-8 -*-

from tramway.inference.base import *
from tramway.inference.smooth_d import smooth_d_neg_posterior
# from TRamWAy.tramway.inference.base import smooth_infer_init
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from collections import OrderedDict


setup = {'name': 'snr',
	#'provides': 'd',
	'arguments': OrderedDict((
		('localization_error',	('-e', dict(type=float, default=0.03, help='localization error'))),
		('diffusivity_prior',	('-d', dict(type=float, default=1., help='prior on the diffusivity'))),
		('jeffreys_prior',	('-j', dict(action='store_true', help="Jeffreys' prior"))),
		('min_diffusivity',	dict(type=float, help='minimum diffusivity value allowed')),
		('max_iter',		dict(type=int, help='maximum number of iterations')))),
	'cell_sampling': 'group'}


def infer_snr(cells, localization_error=0.03, diffusivity_prior=1., jeffreys_prior=None, \
	min_diffusivity=None, max_iter=None, **kwargs):
	# initial values and common calculations
	index, reverse_index, n, dt_mean, D_initial, min_diffusivity, D_bounds = \
		smooth_infer_init(cells, min_diffusivity=min_diffusivity, jeffreys_prior=jeffreys_prior)
	# infer diffusivity D (defined at cells `index`)
	if min_diffusivity is not None:
		kwargs['bounds'] = D_bounds
	if max_iter:
		options = kwargs.get('options', {})
		options['maxiter'] = max_iter
		kwargs['options'] = options
	sle = localization_error * localization_error # sle = squared localization error
	result = minimize(smooth_d_neg_posterior, D_initial, \
		args=(cells, sle, diffusivity_prior, jeffreys_prior, dt_mean, min_diffusivity, reverse_index), \
		**kwargs)
	D = result.x
	# compute diffusivity gradient g (defined at cells `g_index`)
	g_index, g = [], []
	g_defined = np.zeros(len(index), dtype=bool)
	for j, i in enumerate(index):
		gradD = cells.grad(i, D, reverse_index)
		if gradD is not None:
			g_defined[j] = True
			g_index.append(i)
			g.append(gradD[np.newaxis,:])
	g = np.concatenate(g, axis=0)
	# compute mean displacement m and variances V and V_prior (defined at cells `index`)
	sum_pts  = lambda a: np.sum(a, axis=0, keepdims=True)
	sum_dims = lambda a: np.sum(a, axis=1, keepdims=True)
	m, dts, dr, dr2 = [], [], [], []
	for i in index:
		cell = cells[i]
		m.append(np.mean(cell.dr, axis=0, keepdims=True))
		dr.append(sum_pts(cell.dr))
		dr2.append(sum_pts(cell.dr * cell.dr))
		dts.append(cell.dt)
	m   = np.concatenate(m, axis=0)
	dts = np.concatenate(dts)
	n   = n[:,np.newaxis]
	dr  = np.concatenate(dr, axis=0)
	dr2 = np.concatenate(dr2, axis=0)
	V   = sum_dims(dr2 - dr * dr / n) / (n - 1)
	#n_prior   = np.sum(n)    - n
	#dr_prior  = sum_pts(dr)  - dr
	#dr2_prior = sum_pts(dr2) - dr2
	#V_prior   = sum_dims(dr2_prior - dr_prior * dr_prior / n_prior) / (n_prior - 1)
	# compute zeta_total (defined at cells `index`) and zeta_spurious (defined at cells `g_index`)
	sd = np.sqrt(V)
	zeta_total = m / sd
	dt = np.median(dts)
	zeta_spurious = g * dt / sd[g_defined]
	# format the output
	result = pd.DataFrame(
		np.concatenate((n, dr, dr2, D[:,np.newaxis], zeta_total), axis=1),
		index=index,
		columns=['n'] + \
			['dr '+col for col in cells.space_cols] + \
			['dr2 '+col for col in cells.space_cols] + \
			['diffusivity'] + \
			['zeta_total '+col for col in cells.space_cols],
		)
	result = result.join(pd.DataFrame(
		zeta_spurious,
		index=g_index,
		columns=['zeta_spurious '+col for col in cells.space_cols],
		))
	return result

