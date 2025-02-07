import plotly.graph_objects as go

def plot_spectrum(frequencies: list[float], intensities: list[float], x_limits: tuple[float, float] = None, sigma_range: tuple[int, int, int] = None, broaden_type: str = 'lorentzian', resolution: int = 3000) -> go.Figure:
    '''
    Plots a broadened vibrational frequency spectrum with a slider to adjust the broadening parameter (sigma).

    Parameters
    ----------
    intensities : list[float]
        list of IR intensities corresponding to the frequencies.
    frequencies : list
        A list of vibrational frequencies in wave numbers.
    x_limits : tuple[float, float], default=None
        Lowest and highest value shown on the x-axis (frequency)
    sigma_range : tuple[int, int, int], default=None
        Broadening parameters (sigmas) given as a range of integer values (1/cm)
    broaden_type : str, default="lorentzian"
        The type of broadening function to use. Can be "gaussian" or "lorentzian".
    resolution : int, default=3000
        The number of data points in the x-axis for the broadened spectrum.

    Returns
    -------
    go.Figure
        A Plotly figure object representing the broadened vibrational frequency spectrum with a slider.
    '''
