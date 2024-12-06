import plotly.graph_objects as go

def plot_spectrum(frequencies: list[float], intensities: list[float], start: float, end: float, resolution: int = 3000, sigma: float = 20, broaden_type: str = 'lorentzian') -> go.Figure:
    '''
    Plots a broadened vibrational frequency spectrum with a slider to adjust the broadening parameter (sigma).

    Parameters
    ----------
    intensities : list[float]
        list of IR intensities corresponding to the frequencies.
    frequencies : list
        A list of vibrational frequencies in wave numbers.
    start : float
        The starting frequency for the plot.
    end : float
        The ending frequency for the plot.
    resolution : int, default=3000
        The number of data points in the x-axis for the broadened spectrum.
    sigma : float, default=20
        The initial broadening parameter (sigma) for the spectrum.
    broaden_type : str, default="lorentzian"
        The type of broadening function to use. Can be "gaussian" or "lorentzian".

    Returns
    -------
    go.Figure
        A Plotly figure object representing the broadened vibrational frequency spectrum with a slider.
    '''
