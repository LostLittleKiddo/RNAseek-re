# rsa/templatetags/file_tags.py
from django import template
import os
import math

register = template.Library()

@register.filter
def basename(value):
    return os.path.basename(value)

@register.filter
def filter_by_type(queryset, type_name):
    """
    Filter a queryset of ProjectFiles by the 'type' field.
    """
    return queryset.filter(type=type_name)

@register.filter
def filter_by_format(queryset, format_name):
    """
    Filter a queryset of ProjectFiles by the 'file_format' field.
    """
    return queryset.filter(file_format=format_name)

@register.filter
def to_significant_digits(value, digits=4):
    try:
        num = float(value)
        if num == 0:
            return "0.000"
        # Calculate the exponent and adjust for significant digits
        abs_num = abs(num)
        if abs_num == 0:
            return "0.000"
        exponent = math.floor(math.log10(abs_num))
        factor = 10 ** (digits - 1 - exponent)
        rounded = round(num * factor) / factor
        # Format to avoid scientific notation for small numbers
        return f"{rounded:.{digits-1}f}" if exponent < -4 else f"{rounded:.{digits}g}"
    except (ValueError, TypeError):
        return value

@register.filter
def is_number(value):
    try:
        float(value)
        return True
    except (ValueError, TypeError):
        return False