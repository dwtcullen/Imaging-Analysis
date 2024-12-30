"""
Timekeeping toolbox | David Cullen | Jones Lab 2024
"""
import time
def current_time():
    current_time = time.time()  # Use a different variable name
    return current_time

from datetime import datetime, timedelta
def get_eta_string(duration_left: float) -> str:
    eta = datetime.now() + timedelta(seconds=duration_left)
    if eta.date() == datetime.now().date():
        eta_str = "today"
    elif eta.date() == datetime.now().date() + timedelta(days=1):
        eta_str = "tomorrow"
    else:
        eta_str = eta.strftime("%A")
    eta_str += " " + eta.strftime("%I:%M%p").lstrip("0").lower()
    return f"ETA {eta_str}"

def progress(i: int,
             total_groups: int,
             name: tuple[str, str, str],
             start_time,
            ) -> tuple:
    """
    Displays progress information for a loop iteration.

    Parameters:
    i (int): The current iteration number.
    total_groups (int): The total number of groups to iterate over.
    name (Tuple[str, str, str]): A tuple containing metadata information 
                                 about the current iteration (well, row, column).

    Returns:
    progress_report (str): A string containing the time remaining until completion, the
    progress as percentage, the number of total fields, the well ID and field ID.
    """
    current_time = time.time()
    duration_so_far = current_time - start_time
    mean_duration = duration_so_far / i if i > 0 else 0
    duration_left = mean_duration * (total_groups-i)
    
    hours, remainder = divmod(duration_left, 3600)
    minutes, seconds = divmod(remainder, 60)
    duration_left_str = f"{int(hours):02d}h:{int(minutes):02d}m"
    
    current_hours, current_remainder = divmod(duration_so_far, 3600)
    current_minutes, current_seconds = divmod(current_remainder, 60)
    so_far_str = f"{int(current_hours):02d}h:{int(current_minutes):02d}m"
    
    eta_str = get_eta_string(duration_left)
    completion_time = time.strftime('%H:%M:%S', time.localtime(current_time + duration_left))
    percent_progress = ((i)/total_groups)*100
    formatted_progress = "{:.2f}".format(round(percent_progress, 2))
    progress_tuple = (formatted_progress, i, total_groups, name)
    progress_report = f"{eta_str}|{so_far_str}/{duration_left_str}|{progress_tuple[0]}%|{progress_tuple[1]+1}/{progress_tuple[2]}|{name[0]}{name[1]} field {name[2]}"
    return progress_report
