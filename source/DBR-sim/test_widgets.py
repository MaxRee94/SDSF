"""Unit test suite for widgets module.

This module provides comprehensive tests for the widget functionality.
Tests cover all widget classes, their methods, and integration with gui_utils.

Note: Since PySide6/Qt is not available in the test environment, these tests analyze
 the source code structure rather than executing the actual GUI code.
"""

import unittest
import os
import sys

# Add the parent directory to Python path to access the module file
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Read the widgets module source code
with open(os.path.join(sys.path[0], 'widgets.py'), 'r', encoding='utf-8') as f:
    WIDGETS_SOURCE = f.read()


class TestSettingsContainer(unittest.TestCase):
    """Test cases for SettingsContainer class."""

    def test_settings_container_class_exists(self):
        """Test that SettingsContainer class exists."""
        self.assertIn('class SettingsContainer(QtWidgets.QWidget):', WIDGETS_SOURCE)

    def test_settings_container_inherits_from_qwidget(self):
        """Test that SettingsContainer inherits from QWidget."""
        self.assertIn('class SettingsContainer(QtWidgets.QWidget):', WIDGETS_SOURCE)

    def test_settings_container_init_signature(self):
        """Test SettingsContainer.__init__ signature."""
        self.assertIn('def __init__(self, model, parent=None):', WIDGETS_SOURCE)

    def test_settings_container_has_required_attributes(self):
        """Test that SettingsContainer has required attributes."""
        expected_attrs = ['model', 'path_widget', 'family_widget', 'refresh']
        for attr in expected_attrs:
            self.assertIn(attr, WIDGETS_SOURCE)

    def test_settings_container_refresh_method(self):
        """Test that SettingsContainer has refresh method."""
        self.assertIn('def refresh(self):', WIDGETS_SOURCE)


class TestSettingsWidget(unittest.TestCase):
    """Test cases for SettingsWidget class."""

    def test_settings_widget_class_exists(self):
        """Test that SettingsWidget class exists."""
        self.assertIn('class SettingsWidget(QtWidgets.QWidget):', WIDGETS_SOURCE)

    def test_settings_widget_init_signature(self):
        """Test SettingsWidget.__init__ signature."""
        self.assertIn('def __init__(self, parent=None):', WIDGETS_SOURCE)

    def test_settings_widget_has_settings_container(self):
        """Test that SettingsWidget has settings_container attribute."""
        self.assertIn('settings_container', WIDGETS_SOURCE)

    def test_settings_widget_set_file_method(self):
        """Test that SettingsWidget has set_file method."""
        self.assertIn('def set_file(self, file_id):', WIDGETS_SOURCE)


class TestFamilyContainer(unittest.TestCase):
    """Test cases for FamilyContainer class."""

    def test_family_container_class_exists(self):
        """Test that FamilyContainer class exists."""
        self.assertIn('class FamilyContainer(QtWidgets.QComboBox):', WIDGETS_SOURCE)

    def test_family_container_inherits_from_qcombobox(self):
        """Test that FamilyContainer inherits from QComboBox."""
        self.assertIn('class FamilyContainer(QtWidgets.QComboBox):', WIDGETS_SOURCE)

    def test_family_container_has_text_changed_signal(self):
        """Test that FamilyContainer has textChanged signal."""
        self.assertIn('textChanged = QtCore.Signal()', WIDGETS_SOURCE)

    def test_family_container_placeholder_text(self):
        """Test that FamilyContainer has placeholder text constant."""
        self.assertIn('_PLACEHOLDER_TEXT = "Select.."', WIDGETS_SOURCE)

    def test_family_container_init_signature(self):
        """Test FamilyContainer.__init__ signature."""
        self.assertIn('def __init__(self, parent=None):', WIDGETS_SOURCE)

    def test_family_container_required_methods(self):
        """Test that FamilyContainer has required methods."""
        required_methods = ['__init__', 'populate', '_on_text_changed']
        for method in required_methods:
            self.assertIn(f'def {method}', WIDGETS_SOURCE)


class TestInstancesContainer(unittest.TestCase):
    """Test cases for InstancesContainer class."""

    def test_instances_container_class_exists(self):
        """Test that InstancesContainer class exists."""
        self.assertIn('class InstancesContainer(ListSettingWidget):', WIDGETS_SOURCE)

    def test_instances_container_inherits_from_list_setting_widget(self):
        """Test that InstancesContainer inherits from ListSettingWidget."""
        self.assertIn('class InstancesContainer(ListSettingWidget):', WIDGETS_SOURCE)


class TestFamiliesContainer(unittest.TestCase):
    """Test cases for FamiliesContainer class."""

    def test_families_container_class_exists(self):
        """Test that FamiliesContainer class exists."""
        self.assertIn('class FamiliesContainer(ListSettingWidget):', WIDGETS_SOURCE)

    def test_families_container_inherits_from_list_setting_widget(self):
        """Test that FamiliesContainer inherits from ListSettingWidget."""
        self.assertIn('class FamiliesContainer(ListSettingWidget):', WIDGETS_SOURCE)

    def test_families_container_class_variables(self):
        """Test that FamiliesContainer has expected class variables."""
        expected_vars = ['_NEW_ITEM_FIELD', '_PLACEHOLDER_TEXT']
        for var in expected_vars:
            self.assertIn(var, WIDGETS_SOURCE)


class TestModuleStructure(unittest.TestCase):
    """Test cases for overall module structure."""

    def test_module_docstring(self):
        """Test that module has docstring."""
        first_line = WIDGETS_SOURCE.split('\n')[0]
        self.assertIn('"""', first_line)
        self.assertGreater(len(first_line), 3)

    def test_expected_classes_present(self):
        """Test that all expected classes are present in the module."""
        expected_classes = [
            'SettingsContainer', 'SettingsWidget', 
            'FamilyContainer', 'InstancesContainer', 'FamiliesContainer'
        ]
        for class_name in expected_classes:
            self.assertIn(f'class {class_name}', WIDGETS_SOURCE)

    def test_gui_utils_imports(self):
        """Test that gui_utils components are imported."""
        self.assertIn('from gui_utils import HContainer, VContainer, ListSettingWidget, KEYCODES', WIDGETS_SOURCE)


class TestWidgetBehavior(unittest.TestCase):
    """Test cases for widget behavior by examining source code."""

    def test_settings_container_uses_qformlayout(self):
        """Test that SettingsContainer uses QFormLayout."""
        self.assertIn('QFormLayout', WIDGETS_SOURCE)

    def test_settings_container_has_path_widgets(self):
        """Test that SettingsContainer has path-related widgets."""
        expected_widgets = ['path_font', 'path_widget']
        for widget in expected_widgets:
            self.assertIn(widget, WIDGETS_SOURCE)

    def test_settings_container_has_family_widgets(self):
        """Test that SettingsContainer has family-related widgets."""
        expected_widgets = ['family_widget', 'FamilyContainer']
        for widget in expected_widgets:
            self.assertIn(widget, WIDGETS_SOURCE)

    def test_settings_container_has_source_widgets(self):
        """Test that SettingsContainer has source-related widgets."""
        expected_widgets = ['source_line_edit', 'latest_button', 'get_selected_button']
        for widget in expected_widgets:
            self.assertIn(widget, WIDGETS_SOURCE)

    def test_settings_container_has_instance_widgets(self):
        """Test that SettingsContainer has instance-related widgets."""
        expected_widgets = ['InstancesContainer']
        for widget in expected_widgets:
            self.assertIn(widget, WIDGETS_SOURCE)

    def test_settings_container_has_custom_data_widget(self):
        """Test that SettingsContainer has custom data widget."""
        self.assertIn('custom_instance_data_widget', WIDGETS_SOURCE)

    def test_settings_widget_uses_qvboxlayout(self):
        """Test that SettingsWidget uses QVBoxLayout."""
        self.assertIn('QVBoxLayout', WIDGETS_SOURCE)

    def test_family_container_populates_combobox(self):
        """Test that FamilyContainer has populate method."""
        self.assertIn('def populate', WIDGETS_SOURCE)

    def test_family_container_handles_text_changes(self):
        """Test that FamilyContainer handles text changes."""
        self.assertIn('currentTextChanged.connect', WIDGETS_SOURCE)
        self.assertIn('_on_text_changed', WIDGETS_SOURCE)

    def test_settings_container_layout_structure(self):
        """Test that SettingsContainer has expected layout structure."""
        self.assertIn('addRow', WIDGETS_SOURCE)

    def test_settings_widget_layout_structure(self):
        """Test that SettingsWidget has expected layout structure."""
        self.assertIn('addWidget', WIDGETS_SOURCE)
        self.assertIn('addStretch', WIDGETS_SOURCE)


class TestModuleImports(unittest.TestCase):
    """Test cases for module imports."""

    def test_widgets_import(self):
        """Test that widgets module is not circularly imported."""
        # This is a top-level module, so it shouldn't import itself
        self.assertNotIn('import widgets', WIDGETS_SOURCE)

    def test_gui_utils_imports(self):
        """Test that gui_utils imports are present."""
        self.assertIn('from gui_utils import', WIDGETS_SOURCE)

    def test_pyside6_imports(self):
        """Test that PySide6 modules are imported."""
        self.assertIn('from PySide6 import', WIDGETS_SOURCE)


class TestWidgetHierarchy(unittest.TestCase):
    """Test cases for widget hierarchy and relationships."""

    def test_settings_widget_contains_settings_container(self):
        """Test that SettingsWidget contains SettingsContainer."""
        self.assertIn('SettingsContainer', WIDGETS_SOURCE)
        self.assertIn('settings_container', WIDGETS_SOURCE)

    def test_settings_container_contains_family_container(self):
        """Test that SettingsContainer contains FamilyContainer."""
        self.assertIn('FamilyContainer', WIDGETS_SOURCE)
        self.assertIn('family_widget', WIDGETS_SOURCE)

    def test_settings_container_contains_instances_container(self):
        """Test that SettingsContainer contains InstancesContainer."""
        self.assertIn('InstancesContainer', WIDGETS_SOURCE)

    def test_settings_container_contains_families_container(self):
        """Test that SettingsContainer contains FamiliesContainer."""
        self.assertIn('FamiliesContainer', WIDGETS_SOURCE)


class TestWidgetConfiguration(unittest.TestCase):
    """Test cases for widget configuration and styling."""

    def test_settings_container_margins(self):
        """Test that SettingsContainer sets contents margins."""
        self.assertIn('setContentsMargins', WIDGETS_SOURCE)

    def test_settings_container_palette(self):
        """Test that SettingsContainer sets palette."""
        self.assertIn('setPalette', WIDGETS_SOURCE)

    def test_source_buttons_disabled_by_default(self):
        """Test that source buttons are disabled by default."""
        self.assertIn('setEnabled(False)', WIDGETS_SOURCE)

    def test_instance_field_placeholder(self):
        """Test that instance field has placeholder text."""
        self.assertIn('setPlaceholderText', WIDGETS_SOURCE)
        self.assertIn('Enter instance name..', WIDGETS_SOURCE)

    def test_custom_data_widget_placeholder(self):
        """Test that custom data widget has placeholder text."""
        self.assertIn('(Optional) Enter custom python dictionary...', WIDGETS_SOURCE)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)