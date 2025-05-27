# rsa/forms.py
from django import forms
import re

class MultipleFileInput(forms.ClearableFileInput):
    allow_multiple_selected = True

class MultipleFileField(forms.FileField):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("widget", MultipleFileInput())
        super().__init__(*args, **kwargs)

    def clean(self, data, initial=None):
        single_file_clean = super().clean
        if isinstance(data, (list, tuple)):
            result = [single_file_clean(d, initial) for d in data]
        else:
            result = [single_file_clean(data, initial)]
        return result

class RNAseekForm(forms.Form):
    project_name = forms.CharField(
        max_length=255,
        required=True,
        label="Project Name",
        widget=forms.TextInput(attrs={
            'placeholder': 'Enter a name for your analysis',
            'class': 'mt-1 block w-full px-4 py-3 border border-gray-300 rounded-lg shadow-sm focus:ring-emerald-500 focus:border-emerald-500 sm:text-sm transition-all duration-300'
        })
    )
    genome_of_interest = forms.ChoiceField(
        choices=[
            ('', 'Select a genome'),
            ('yeast', 'Yeast (Saccharomyces cerevisiae, R64-1-1)'),
            ('human', 'Human (Homo sapiens, GRCh38)'),
            ('mouse', 'Mouse (Mus musculus, GRCm39)')
        ],
        required=True,
        label="Genome of Interest",
        widget=forms.Select(attrs={
            'class': 'mt-1 block w-full px-4 py-3 border border-gray-300 rounded-lg shadow-sm focus:ring-emerald-500 focus:border-emerald-500 sm:text-sm transition-all duration-300'
        })
    )
    sequencing_type = forms.ChoiceField(
        choices=[
            ('single', 'Single-End'),
            ('paired', 'Paired-End')
        ],
        required=True,
        label="Sequencing Type",
        initial='single',
        widget=forms.RadioSelect(attrs={
            'class': 'h-4 w-4 text-emerald-600 border-gray-300 focus:ring-emerald-500'
        })
    )
    pvalue_cutoff = forms.FloatField(
        required=False,
        min_value=0.0,
        max_value=1.0,
        initial=0.05,
        label="P-value Cutoff",
        widget=forms.NumberInput(attrs={
            'step': '0.01',
            'placeholder': 'Enter p-value cutoff (0 to 1)',
            'class': 'mt-1 block w-full px-4 py-3 border border-gray-300 rounded-lg shadow-sm focus:ring-emerald-500 focus:border-emerald-500 sm:text-sm transition-all duration-300'
        })
    )
    files = MultipleFileField(
        required=True,
        label="Upload Files",
        widget=MultipleFileInput(attrs={
            'multiple': True,
            'accept': '.fastq,.fastq.gz',
            'class': 'hidden',
            'id': 'file-input'
        })
    )

    def clean_files(self):
        files = self.cleaned_data.get('files', [])
        allowed_extensions = ['fastq', 'fq', 'gz']
        max_size = 70 * 1024 * 1024 * 1024  # 70GB
        sequencing_type = self.cleaned_data.get('sequencing_type')

        if sequencing_type == 'single' and len(files) < 2:
            raise forms.ValidationError("For Single-End sequencing, you must upload at least 2 files.")
        if sequencing_type == 'paired' and len(files) < 4:
            raise forms.ValidationError("For Paired-End sequencing, you must upload at least 4 files.")

        for file in files:
            ext = file.name.split('.')[-1].lower()
            if ext not in allowed_extensions:
                raise forms.ValidationError(f"File {file.name} has an invalid extension. Allowed: .fastq, .fastq.gz")
            if file.size > max_size:
                raise forms.ValidationError(f"File {file.name} is too large. Maximum size is 70GB.")

        if sequencing_type == 'paired':
            if len(files) % 2 != 0:
                raise forms.ValidationError("For Paired-End sequencing, you must upload an even number of files (e.g., pairs of forward and reverse reads).")
            
            # Check for paired files using a more robust regex
            sample_names = {}
            for file in files:
                # Match sample name and direction (R1 or R2, case-insensitive, allowing underscores, dots, or numbers)
                match = re.match(r'^(.*?)(?:_R[12])(?:_|\.|\d+|$)', file.name, re.IGNORECASE)
                if not match:
                    raise forms.ValidationError(f"File {file.name} does not contain 'R1' or 'R2' (case-insensitive) in the filename for paired-end sequencing.")
                sample_name = match.group(1)
                direction = 'r1' if 'R1' in file.name.upper() else 'r2'
                if direction not in ['r1', 'r2']:
                    raise forms.ValidationError(f"File {file.name} must indicate R1 or R2 in the filename for paired-end sequencing.")
                
                if sample_name not in sample_names:
                    sample_names[sample_name] = {'r1': None, 'r2': None}
                sample_names[sample_name][direction] = file.name

            # Validate that each sample has both R1 and R2
            for sample_name, directions in sample_names.items():
                if not (directions['r1'] and directions['r2']):
                    raise forms.ValidationError(f"Sample {sample_name} is missing a matching pair. Each sample must have both R1 and R2 files.")

        return files

    def clean(self):
        cleaned_data = super().clean()
        return cleaned_data

class DeseqMetadataForm(forms.Form):
    condition1 = forms.CharField(
        max_length=50,
        required=True,
        label="Condition 1",
        initial="control",
        widget=forms.TextInput(attrs={
            'placeholder': 'e.g., control',
            'class': 'mt-1 block w-full px-4 py-3 border border-gray-300 rounded-lg shadow-sm focus:ring-emerald-500 focus:border-emerald-500 sm:text-sm transition-all duration-300'
        })
    )
    condition2 = forms.CharField(
        max_length=50,
        required=True,
        label="Condition 2",
        initial="treatment",
        widget=forms.TextInput(attrs={
            'placeholder': 'e.g., treatment',
            'class': 'mt-1 block w-full px-4 py-3 border border-gray-300 rounded-lg shadow-sm focus:ring-emerald-500 focus:border-emerald-500 sm:text-sm transition-all duration-300'
        })
    )

    def __init__(self, *args, **kwargs):
        files = kwargs.pop('files', [])
        sequencing_type = kwargs.pop('sequencing_type', 'single')
        super().__init__(*args, **kwargs)

        # Dynamically add fields for each file or sample
        self.sample_names = []
        if sequencing_type == 'single':
            for file in files:
                name = file.name.split('.fastq')[0]
                self.sample_names.append(name)
                self.fields[f'condition_{name}'] = forms.ChoiceField(
                    choices=[('', 'Select condition'), ('condition1', 'Condition 1'), ('condition2', 'Condition 2')],
                    required=True,
                    label=f"Condition for {name}",
                    widget=forms.Select(attrs={
                        'class': 'mt-1 block w-full px-4 py-3 border border-gray-300 rounded-lg shadow-sm focus:ring-emerald-500 focus:border-emerald-500 sm:text-sm transition-all duration-300'
                    })
                )
        else:  # paired
            sample_names = {}
            for file in files:
                # More robust regex for paired-end sample names
                match = re.match(r'^(.*?)(?:_R[12])(?:_|\.|\d+|$)', file.name, re.IGNORECASE)
                if match:
                    sample_name = match.group(1)
                    if sample_name not in sample_names:
                        sample_names[sample_name] = True
                        self.sample_names.append(sample_name)
                        self.fields[f'condition_{sample_name}'] = forms.ChoiceField(
                            choices=[('', 'Select condition'), ('condition1', 'Condition 1'), ('condition2', 'Condition 2')],
                            required=True,
                            label=f"Condition for {sample_name}",
                            widget=forms.Select(attrs={
                                'class': 'mt-1 block w-full px-4 py-3 border border-gray-300 rounded-lg shadow-sm focus:ring-emerald-500 focus:border-emerald-500 sm:text-sm transition-all duration-300'
                            })
                        )

    def clean(self):
        cleaned_data = super().clean()
        condition1 = cleaned_data.get('condition1')
        condition2 = cleaned_data.get('condition2')

        if condition1 and condition2 and condition1 == condition2:
            raise forms.ValidationError("Condition 1 and Condition 2 must be different.")

        # Ensure all sample conditions are selected
        for field_name, field in self.fields.items():
            if field_name.startswith('condition_') and not cleaned_data.get(field_name):
                raise forms.ValidationError(f"Please select a condition for {field.label}.")

        return cleaned_data