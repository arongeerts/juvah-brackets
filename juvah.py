# Load the Python Standard and DesignScript Libraries
import math
import sys
# noinspection PyUnresolvedReferences
import clr

clr.AddReference('ProtoGeometry')
clr.AddReference("RevitAPI")

# noinspection PyUnresolvedReferences
from Autodesk.Revit.DB import *
from Autodesk.DesignScript.Geometry import Curve, Point

# The inputs to this node will be stored as a list in the IN variables.
# noinspection PyUnresolvedReferences
dataEnteringNode = IN

# Place your code below this line
ducts = dataEnteringNode[0]

BRACKET_FAMILY = "VE_Bracket_MEPcontent_CAIROX_OBMC"
CHANNEL_FAMILY = "JUV_57.0_DUA_UN_Beugel rechthoekig kanaal"
FEET_TO_MM = 304.8

def get_parameters(duct):
    parameters = {}
    for param in duct.Parameters:
        storage_type = param.StorageType

        # Get the parameter value based on its storage type
        if storage_type == StorageType.String:
            param_value = param.AsString()  # For string parameters
        elif storage_type == StorageType.Integer:
            param_value = param.AsInteger()  # For integer parameters
        elif storage_type == StorageType.Double:
            param_value = param.AsDouble()  # For double parameters (often stored in feet)
        elif storage_type == StorageType.ElementId:
            param_value = param.AsElementId().IntegerValue  # For ElementId parameters (e.g., materials)
        else:
            param_value = "Unknown Type"
        parameters[param.Definition.Name] = param_value
    return parameters


def get_diameter(duct):
    parameters = get_parameters(duct)
    return float(parameters["D1 Description"])


def is_round(duct):
    for param in duct.Parameters:
        if param.Definition.Name == "Diameter":
            return True
    return False


def get_duct_orientation(duct):
    curve = duct.Location.Curve

    # Get the start and end points of the duct
    start_point = curve.GetEndPoint(0)  # Starting point
    end_point = curve.GetEndPoint(1)  # Ending point

    # Get the XY coordinates for angle calculation
    x = end_point.X - start_point.X
    y = end_point.Y - start_point.Y

    # Calculate the angle in degrees
    angle_xy = math.degrees(math.atan2(y, x)) % 360

    # Diagonal ducts get a 90 degree offset, reason unknown
    result = round(angle_xy)
    if result not in [0, 90, 180, 270, 360]:
        result += 90

    return result


def get_length(duct):
    return duct.Location.Curve.Length * FEET_TO_MM


def point_from_coordinate(coordinate):
    # Convert from feet to mm
    return Point.ByCoordinates(
        coordinate.X * FEET_TO_MM,
        coordinate.Y * FEET_TO_MM,
        coordinate.Z * FEET_TO_MM
    )


def get_points_on_duct(duct, curve, max_spacing, offset):
    length = get_length(duct)
    span = length - (2 * offset)
    amount = math.ceil(span / max_spacing) + 1
    spacing = span / (amount - 1)
    points = []

    for value in range(int(amount)):
        ratio = (offset + (value * spacing)) / length
        coordinates = duct.Location.Curve.Evaluate(ratio, True)
        points.append(point_from_coordinate(coordinates))
    return points


def handle_diameter_less_than_200(duct):
    diameter = get_diameter(duct)
    if get_length(duct) <= 1000:
        coordinates = duct.Location.Curve.Evaluate(0.5, True)
        points = [point_from_coordinate(coordinates)]
    else:
        points = get_points_on_duct(duct, duct.Location, 1500, 500)

    return [[BRACKET_FAMILY, get_bracket_family_type(diameter)], [points, get_duct_orientation(duct)]]


def handle_diameter_between_200_and_450(duct):
    diameter = get_diameter(duct)
    if get_length(duct) <= 1000:
        coordinates = duct.Location.Curve.Evaluate(0.5, True)
        points = [point_from_coordinate(coordinates)]
    else:
        points = get_points_on_duct(duct, duct.Location, 1000, 500)

    return [[BRACKET_FAMILY, get_bracket_family_type(diameter)], [points, get_duct_orientation(duct)]]


def get_bracket_family_type(diameter):
    t = str(int(diameter))
    if t not in ["125", "160", "200", "250"]:
        return "Default"
    return t


def handle_diameter_greater_than_450(duct):
    # Access the width parameter of the rectangular duct (typically called "Width" or "Diameter Width")
    width = get_diameter(duct)

    return get_canals(duct, width)


def handle_rectangular(duct):
    width_param = duct.LookupParameter("Width")
    width = width_param.AsDouble() * FEET_TO_MM

    return get_canals(duct, width)


def get_canals(duct, width):
    curve = duct.Location.Curve

    if get_length(duct) <= 1500:
        coordinates = duct.Location.Curve.Evaluate(0.5, True)
        points = [point_from_coordinate(coordinates)]
        return [
            points,
            width,
            get_duct_orientation(duct) + 90 % 360
        ]
    else:
        points = get_points_on_duct(duct, curve, 1500, 500)

        return [
            points,
            width,
            get_duct_orientation(duct) + 90 % 360
        ]


def is_vertical(duct):
    curve = duct.Location.Curve

    # Get the start and end points of the duct
    start_point = curve.GetEndPoint(0)  # Starting point
    end_point = curve.GetEndPoint(1)  # Ending point

    return abs(start_point.X - end_point.X) < 0.0001 and abs(start_point.Y - end_point.Y) < 0.0001


def main():
    brackets = []
    canals = []
    for duct in ducts:
        duct = UnwrapElement(duct)
        if is_vertical(duct):
            continue
        if is_round(duct):
            diameter = get_diameter(duct)
            if not diameter:
                continue
            if diameter <= 200:
                brackets.append(handle_diameter_less_than_200(duct))
            elif diameter > 200 and diameter <= 450:
                brackets.append(handle_diameter_between_200_and_450(duct))
            else:
                canals.append(handle_diameter_greater_than_450(duct) or [])
        else:
            canals.append(handle_rectangular(duct) or [])

        brackets = [b for b in brackets if b is not None]
        canals = [c for c in canals if c is not None]

    return brackets, canals


OUT = main()
