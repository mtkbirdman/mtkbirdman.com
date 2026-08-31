import numpy as np

RHO0 = 1.225
G = 9.80665


def atmosphere(case, altitude_m):
    """Return density and ground-axis wind at one altitude or an altitude array."""
    settings = case["atmosphere"]
    altitude = np.maximum(0.0, np.asarray(altitude_m, dtype=float))

    if settings["model"] == "constant":
        rho = np.full_like(altitude, float(settings.get("density_kg_m3", RHO0)))
    elif settings["model"] == "isa":
        # Troposphere only; more than sufficient for the intended winch-launch range.
        temperature = 288.15 - 0.0065 * altitude
        pressure = 101325.0 * (temperature / 288.15) ** 5.2558797
        rho = pressure / (287.05287 * temperature)
    else:
        raise ValueError(f"Unsupported atmosphere model: {settings['model']}")

    wind_x = float(settings.get("wind_x_m_s", 0.0)) + float(
        settings.get("wind_x_gradient_per_m", 0.0)
    ) * altitude
    wind_z = float(settings.get("wind_z_m_s", 0.0)) + float(
        settings.get("wind_z_gradient_per_m", 0.0)
    ) * altitude

    if altitude.ndim == 0:
        return float(rho), float(wind_x), float(wind_z)
    return rho, wind_x, wind_z
