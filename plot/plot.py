import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame

def plot_experiments_results(data: DataFrame, type_data: str, save_fig: bool = False, fig_name: str = None) -> None:
    """Creates charts with input DataFrame.

    Args:
        data (DataFrame): The data obtained in experiments.
        type (str): If we want to plot cnot infidelity experiment results or qec infidelity experiment results.
        save_fig (bool, optional): If we want to save the plot. Defaults to False.
        fig_name (str, optional): If we want to add a specific name to the plot. Defaults to None.

    Raises:
        TypeError: If data is not equal a DataFrame.
        TypeError: If fig_name is not a str.
    """
    if isinstance(data, DataFrame):
        
        if type_data == "infidelity":
            
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
        elif type_data == "qec":
            
            x1, y1 = np.polyfit(x = data["Process infidelity"], y = data["p_logical_error"], deg = 1)
            x2, y2 = np.polyfit(x = data["Average gate infidelity"], y = data["p_logical_error"], deg = 1)
            plt.style.use('ggplot')
            fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2)
            fig.set_size_inches([12, 9])
            ax1.plot(data["Process infidelity"], x1 * data["Process infidelity"] + y1, color = "blue", label = "Linear regression")
            ax1.scatter(data["Process infidelity"], data["p_logical_error"], marker = "o", color = "black", label = "Data points")
            ax1.set_xlabel("Process infidelity")
            ax1.set_ylabel("Logical error probability")
            ax1.legend()
            ax2.plot(data["Average gate infidelity"], x2 * data["Average gate infidelity"] + y2, color = "blue", label = "Linear regression")
            ax2.scatter(data["Average gate infidelity"], data["p_logical_error"], marker = "o", color = "black", label = "Data points")
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
        else:
            
            raise ValueError("The type_data must be equal to infidelity or qec!")
    else:
        
        raise TypeError("The input is not a DataFrame!")

# def plot_code_experiments_results(data: DataFrame, save_fig: bool = False, fig_name: str = None) -> None:
#     """Creates charts with input DataFrame.

#     Args:
#         data1 (DataFrame): The data from experiments.
#         save_fig (bool, optional): If we want to save the figure. Defaults to False.
#         fig_name (str, optional): If we want to add a specific name to the plot. Defaults to None.

#     Raises:
#         TypeError: If data is not equal to a DataFrame.
#         TypeError: If fig_name is not a str.
#     """
#     if isinstance(data, DataFrame):
        
#         x1, y1 = np.polyfit(x = data["Process infidelity"], y = data["p_logical_error"], deg = 1)
#         x2, y2 = np.polyfit(x = data["Average gate infidelity"], y = data["p_logical_error"], deg = 1)
#         plt.style.use('ggplot')
#         fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2)
#         fig.set_size_inches([12, 9])
#         ax1.plot(data["Process infidelity"], x1 * data["Process infidelity"] + y1, color = "blue", label = "Linear regression")
#         ax1.scatter(data["Process infidelity"], data["p_logical_error"], marker = "o", color = "black", label = "Data points")
#         ax1.set_xlabel("Process infidelity")
#         ax1.set_ylabel("Logical error probability")
#         ax1.legend()
#         ax2.plot(data["Average gate infidelity"], x2 * data["Average gate infidelity"] + y2, color = "blue", label = "Linear regression")
#         ax2.scatter(data["Average gate infidelity"], data["p_logical_error"], marker = "o", color = "black", label = "Data points")
#         ax2.set_xlabel("Average gate infidelity")
#         ax2.legend()
#         if save_fig:
        
#             if fig_name is None:
            
#                 fig.savefig("figures/plot.jpg", quality = 100)
#             else:
                
#                 if isinstance(fig_name, str):
                    
#                     fig.savefig(f"figures/{fig_name}.jpg", quality = 100)
#                 else:
                    
#                     raise TypeError("The input is not a str!")
        
#     else:
        
#         raise TypeError("The input is not equal to a DataFrame!")