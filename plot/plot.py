import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame

def plot_infidelity_experiments_results(data: DataFrame, save_fig: bool = False, fig_name: str = None) -> None:
    """Creates charts with input DataFrame.

    Args:
        data (DataFrame): The data obtained in experiments.
        save_fig (bool, optional): If we want to save the plot. Defaults to False.
        fig_name (str, optional): If we want to add a specific name to the plot. Defaults to None.

    Raises:
        TypeError: If data is not equal a DataFrame.
        TypeError: If fig_name is not a str.
    """
    if isinstance(data, DataFrame):
        
        x1, y1 = np.polyfit(x = data["Process infidelity"], y = data["p_not_00"], deg = 1)
        x2, y2 = np.polyfit(x = data["Average gate infidelity"], y = data["p_not_00"], deg = 1)
        plt.style.use('ggplot')
        fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2)
        fig.set_size_inches([12, 9])
        ax1.plot(data["Process infidelity"], x1 * data["Process infidelity"] + y1, color = "blue", label = "Linear regression")
        ax1.scatter(data["Process infidelity"], data["p_not_00"], marker = "o", color = "black", label = "Data points")
        ax1.set_xlabel("Process infidelity")
        ax1.set_ylabel("Probability output is not 00")
        ax1.legend()
        ax2.plot(data["Average gate infidelity"], x2 * data["Average gate infidelity"] + y2, color = "blue", label = "Linear regression")
        ax2.scatter(data["Average gate infidelity"], data["p_not_00"], marker = "o", color = "black", label = "Data points")
        ax2.set_xlabel("Average gate infidelity")
        ax2.legend()
        if save_fig:
            
            if fig_name is None:
                
                fig.savefig("figures/plot.jpg", quality = 100)
            else:
                
                if isinstance(fig_name, str):
                    
                    fig.savefig(f"figures/{fig_name}.jpg", quality = 100)
                else:
                    
                    raise TypeError("The input is not a str!")
        plt.show()
    else:
        
        raise TypeError("The input is not a DataFrame!")

def plot_code_experiments_results(data1: DataFrame, data2: DataFrame, save_fig: bool = False, fig_name: str = None) -> None:
    """Creates charts with input DataFrame.

    Args:
        data1 (DataFrame): The data from infidelity experiments.
        data2 (DataFrame): The data from code and logical errors experiments.
        save_fig (bool, optional): If we want to save the figure. Defaults to False.
        fig_name (str, optional): If we want to add a specific name to the plot. Defaults to None.

    Raises:
        TypeError: If data1 and data2 are not equal to a DataFrame.
        TypeError: If fig_name is not a str.
    """
    if isinstance(data1, DataFrame) and isinstance(data2, DataFrame):
        
        x1, y1 = np.polyfit(x = data1["Process infidelity"], y = data2["p_logical_error"], deg = 1)
        x2, y2 = np.polyfit(x = data1["Average gate infidelity"], y = data2["p_logical_error"], deg = 1)
        plt.style.use('ggplot')
        fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2)
        fig.set_size_inches([12, 9])
        ax1.plot(data1["Process infidelity"], x1 * data1["Process infidelity"] + y1, color = "blue", label = "Linear regression")
        ax1.scatter(data1["Process indidelity"], data2["p_logical_error"], marker = "o", color = "black", label = "Data points")
        ax1.set_xlabel("Process infidelity")
        ax1.set_ylabel("Logical error probability")
        ax1.legend()
        ax2.plot(data1["Average gate infidelity"], x2 * data1["Average gate infidelity"] + y2, color = "blue", label = "Linear regression")
        ax2.scatter(data1["Average gate indidelity"], data2["p_logical_error"], marker = "o", color = "black", label = "Data points")
        ax2.set_xlabel("Average gate infidelity")
        ax2.legend()
        if save_fig:
        
            if fig_name is None:
            
                fig.savefig("figures/plot.jpg", quality = 100)
            else:
                
                if isinstance(fig_name, str):
                    
                    fig.savefig(f"figures/{fig_name}.jpg", quality = 100)
                else:
                    
                    raise TypeError("The input is not a str!")
        plt.show()
        
    else:
        
        raise TypeError("The inputs are not equal to DataFrame!")