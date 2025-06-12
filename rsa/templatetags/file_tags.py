# rsa/templatetags/file_tags.py
from django import template
import os

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
        # Convert to float and handle scientific notation
        num = float(value)
        if num == 0:
            return "0.000"
        # Use format specifier for significant digits
        return f"{num:.{digits}g}"
    except (ValueError, TypeError):
        return value  # Return original value if not a number
    
@register.filter
def is_number(value):
    try:
        float(value)
        return True
    except (ValueError, TypeError):
        return False