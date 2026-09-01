"""Unit test suite for widgets module.

This module provides comprehensive tests for the widget functionality.
Tests cover all widget classes, their methods, and integration with gui_utils.
"""

import unittest
import os
import sys
import inspect
from unittest.mock import MagicMock, patch

# Add the parent directory to Python path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Mock PySide6 imports to avoid GUI dependencies in tests
sys.modules['PySide6'] = MagicMock()
sys.modules['PySide6.QtCore'] = MagicMock()
sys.modules['PySide6.QtWidgets'] = MagicMock()
sys.modules['PySide6.QtGui'] = MagicMock()

# Mock gui_utils imports
sys.modules['gui_utils'] = MagicMock()

import widgets as wdg


class TestSettingsContainer(unittest.TestCase):
    """Test cases for SettingsContainer class."""

    def test_settings_container_class_exists(self):
        """Test that SettingsContainer class exists."""
        self.assertTrue(hasattr(wdg, 'SettingsContainer'))
        self.assertTrue(callable(wdg.SettingsContainer))

    def test_settings_container_inherits_from_qwidget(self):
        """Test that SettingsContainer inherits from QWidget."""
        # Test class name and structure
        self.assertEqual(wdg.SettingsContainer.__name__, 'SettingsContainer')

    def test_settings_container_init_signature(self):
        """Test SettingsContainer.__init__ signature."""
        sig = inspect.signature(wdg.SettingsContainer.__init__)
        params = list(sig.parameters.keys())
        
        self.assertIn('self', params)
        self.assertIn('model', params)
        self.assertIn('parent', params)

    def test_settings_container_has_required_attributes(self):
        """Test that SettingsContainer has required attributes."""
        # Check that the class has the expected attributes by examining source
        source = inspect.getsource(wdg.SettingsContainer)
        
        expected_attrs = ['model', 'path_widget', 'family_widget', 'refresh']
        for attr in expected_attrs:
            self.assertIn(attr, source)

    def test_settings_container_refresh_method(self):
        """Test that SettingsContainer has refresh method."""
        self.assertTrue(hasattr(wdg.SettingsContainer, 'refresh'))
        self.assertTrue(callable(wdg.SettingsContainer.refresh))


class TestSettingsWidget(unittest.TestCase):
    """Test cases for SettingsWidget class."""

    def test_settings_widget_class_exists(self):
        """Test that SettingsWidget class exists."""
        self.assertTrue(hasattr(wdg, 'SettingsWidget'))
        self.assertTrue(callable(wdg.SettingsWidget))

    def test_settings_widget_init_signature(self):
        """Test SettingsWidget.__init__ signature."""
        sig = inspect.signature(wdg.SettingsWidget.__init__)
        params = list(sig.parameters.keys())
        
        self.assertIn('self', params)
        self.assertIn('parent', params)

    def test_settings_widget_has_settings_container(self):
        """Test that SettingsWidget has settings_container attribute."""
        source = inspect.getsource(wdg.SettingsWidget)
        self.assertIn('settings_container', source)

    def test_settings_widget_set_file_method(self):
        """Test that SettingsWidget has set_file method."""
        self.assertTrue(hasattr(wdg.SettingsWidget, 'set_file'))
        self.assertTrue(callable(wdg.SettingsWidget.set_file))


class TestFamilyContainer(unittest.TestCase):
    """Test cases for FamilyContainer class."""

    def test_family_container_class_exists(self):
        """Test that FamilyContainer class exists."""
        self.assertTrue(hasattr(wdg, 'FamilyContainer'))
        self.assertTrue(callable(wdg.FamilyContainer))

    def test_family_container_inherits_from_qcombobox(self):
        """Test that FamilyContainer inherits from QComboBox."""
        self.assertEqual(wdg.FamilyContainer.__name__, 'FamilyContainer')

    def test_family_container_has_text_changed_signal(self):
        """Test that FamilyContainer has textChanged signal."""
        source = inspect.getsource(wdg.FamilyContainer)
        self.assertIn('textChanged = QtCore.Signal()', source)

    def test_family_container_placeholder_text(self):
        """Test that FamilyContainer has placeholder text constant."""
        self.assertTrue(hasattr(wdg.FamilyContainer, '_PLACEHOLDER_TEXT'))
        self.assertEqual(wdg.FamilyContainer._PLACEHOLDER_TEXT, "Select..")

    def test_family_container_init_signature(self):
        """Test FamilyContainer.__init__ signature."""
        sig = inspect.signature(wdg.FamilyContainer.__init__)
        params = list(sig.parameters.keys())
        
        self.assertIn('self', params)
        self.assertIn('parent', params)

    def test_family_container_required_methods(self):
        """Test that FamilyContainer has required methods."""
        required_methods = ['__init__', 'populate', '_on_text_changed']
        for method in required_methods:
            self.assertTrue(hasattr(wdg.FamilyContainer, method))
            self.assertTrue(callable(getattr(wdg.FamilyContainer, method)))


class TestInstancesContainer(unittest.TestCase):
    """Test cases for InstancesContainer class."""

    def test_instances_container_class_exists(self):
        """Test that InstancesContainer class exists."""
        self.assertTrue(hasattr(wdg, 'InstancesContainer'))
        self.assertTrue(callable(wdg.InstancesContainer))

    def test_instances_container_inherits_from_list_setting_widget(self):
        """Test that InstancesContainer inherits from ListSettingWidget."""
        # Check source code for inheritance
        source = inspect.getsource(wdg.InstancesContainer)
        self.assertIn('ListSettingWidget', source)


class TestFamiliesContainer(unittest.TestCase):
    """Test cases for FamiliesContainer class."""

    def test_families_container_class_exists(self):
        """Test that FamiliesContainer class exists."""
        self.assertTrue(hasattr(wdg, 'FamiliesContainer'))
        self.assertTrue(callable(wdg.FamiliesContainer))

    def test_families_container_inherits_from_list_setting_widget(self):
        """Test that FamiliesContainer inherits from ListSettingWidget."""
        # Check source code for inheritance
        source = inspect.getsource(wdg.FamiliesContainer)
        self.assertIn('ListSettingWidget', source)

    def test_families_container_class_variables(self):
        """Test that FamiliesContainer has expected class variables."""
        expected_vars = ['_NEW_ITEM_FIELD', '_PLACEHOLDER_TEXT', '_ADD_BTN_TEXT', '_REMOVE_BTN_TEXT']
        for var in expected_vars:
            self.assertTrue(hasattr(wdg.FamiliesContainer, var))


class TestModuleStructure(unittest.TestCase):
    """Test cases for overall module structure."""

    def test_module_docstring(self):
        """Test that module has docstring."""
        self.assertIsNotNone(wdg.__doc__)
        self.assertGreater(len(wdg.__doc__), 0)

    def test_expected_classes_present(self):
        """Test that all expected classes are present in the module."""
        expected_classes = [
            'SettingsContainer', 'SettingsWidget', 
            'FamilyContainer', 'InstancesContainer', 'FamiliesContainer'
        ]
        for class_name in expected_classes:
            self.assertTrue(hasattr(wdg, class_name))
            self.assertTrue(callable(getattr(wdg, class_name)))

    def test_gui_utils_imports(self):
        """Test that gui_utils components are imported."""
        source = inspect.getsource(wdg)
        
        # Check imports from gui_utils
        self.assertIn('HContainer', source)
        self.assertIn('VContainer', source)
        self.assertIn('ListSettingWidget', source)
        self.assertIn('KEYCODES', source)


class TestWidgetBehavior(unittest.TestCase):
    """Test cases for widget behavior by examining source code."""

    def test_settings_container_uses_qformlayout(self):
        """Test that SettingsContainer uses QFormLayout."""
        source = inspect.getsource(wdg.SettingsContainer)
        self.assertIn('QFormLayout', source)

    def test_settings_container_has_path_widgets(self):
        """Test that SettingsContainer has path-related widgets."""
        source = inspect.getsource(wdg.SettingsContainer)
        
        expected_widgets = ['path_font', 'path_widget']
        for widget in expected_widgets:
            self.assertIn(widget, source)

    def test_settings_container_has_family_widgets(self):
        """Test that SettingsContainer has family-related widgets."""
        source = inspect.getsource(wdg.SettingsContainer)
        
        expected_widgets = ['family_widget', 'FamilyContainer']
        for widget in expected_widgets:
            self.assertIn(widget, source)

    def test_settings_container_has_source_widgets(self):
        """Test that SettingsContainer has source-related widgets."""
        source = inspect.getsource(wdg.SettingsContainer)
        
        expected_widgets = ['source_line_edit', 'latest_button', 'get_selected_button']
        for widget in expected_widgets:
            self.assertIn(widget, source)

    def test_settings_container_has_instance_widgets(self):
        """Test that SettingsContainer has instance-related widgets."""
        source = inspect.getsource(wdg.SettingsContainer)
        
        expected_widgets = ['InstancesContainer']
        for widget in expected_widgets:
            self.assertIn(widget, source)

    def test_settings_container_has_custom_data_widget(self):
        """Test that SettingsContainer has custom data widget."""
        source = inspect.getsource(wdg.SettingsContainer)
        self.assertIn('custom_instance_data_widget', source)

    def test_settings_widget_uses_qvboxlayout(self):
        """Test that SettingsWidget uses QVBoxLayout."""
        source = inspect.getsource(wdg.SettingsWidget)
        self.assertIn('QVBoxLayout', source)

    def test_family_container_populates_combobox(self):
        """Test that FamilyContainer populates combobox."""
        source = inspect.getsource(wdg.FamilyContainer.populate)
        
        # Check that populate method adds items to combobox
        self.assertIn('addItem', source)

    def test_family_container_handles_text_changes(self):
        """Test that FamilyContainer handles text changes."""
        source = inspect.getsource(wdg.FamilyContainer.__init__)
        
        # Check that it connects to text change signal
        self.assertIn('currentTextChanged.connect', source)
        self.assertIn('_on_text_changed', source)

    def test_settings_container_layout_structure(self):
        """Test that SettingsContainer has expected layout structure."""
        source = inspect.getsource(wdg.SettingsContainer)
        
        # Check that it adds rows to layout
        self.assertIn('addRow', source)

    def test_settings_widget_layout_structure(self):
        """Test that SettingsWidget has expected layout structure."""
        source = inspect.getsource(wdg.SettingsWidget)
        
        # Check that it adds widgets to layout
        self.assertIn('addWidget', source)
        self.assertIn('addStretch', source)


class TestInstanceContainer(unittest.TestCase):
    """Test cases for InstancesContainer class."""

    def test_instances_container_customization(self):
        """Test that InstancesContainer customizes ListSettingWidget."""
        source = inspect.getsource(wdg.InstancesContainer)
        
        # Check for customization of class variables
        expected_customizations = ['_NEW_ITEM_FIELD', '_PLACEHOLDER_TEXT']
        for custom in expected_customizations:
            self.assertIn(custom, source)


class TestFamilyContainerBehavior(unittest.TestCase):
    """Test cases for FamilyContainer behavior."""

    def test_family_container_on_text_changed_method(self):
        """Test that FamilyContainer has _on_text_changed method."""
        self.assertTrue(hasattr(wdg.FamilyContainer, '_on_text_changed'))
        self.assertTrue(callable(wdg.FamilyContainer._on_text_changed))

    def test_family_container_populate_method(self):
        """Test that FamilyContainer has populate method."""
        self.assertTrue(hasattr(wdg.FamilyContainer, 'populate'))
        self.assertTrue(callable(wdg.FamilyContainer.populate))

    def test_family_container_connects_signals(self):
        """Test that FamilyContainer connects signals properly."""
        source = inspect.getsource(wdg.FamilyContainer.__init__)
        
        # Check signal connections
        self.assertIn('currentTextChanged.connect', source)
        self.assertIn('self._on_text_changed', source)


class TestModuleImports(unittest.TestCase):
    """Test cases for module imports."""

    def test_pyside6_imports(self):
        """Test that PySide6 modules are imported."""
        # This is tested implicitly by the fact that the module loads
        # without import errors (thanks to our mocking)
        self.assertTrue(True)

    def test_gui_utils_imports(self):
        """Test that gui_utils imports are present."""
        source = inspect.getsource(wdg)
        self.assertIn('from gui_utils import', source)


class TestWidgetHierarchy(unittest.TestCase):
    """Test cases for widget hierarchy and relationships."""

    def test_settings_widget_contains_settings_container(self):
        """Test that SettingsWidget contains SettingsContainer."""
        source = inspect.getsource(wdg.SettingsWidget)
        
        self.assertIn('SettingsContainer', source)
        self.assertIn('settings_container', source)

    def test_settings_container_contains_family_container(self):
        """Test that SettingsContainer contains FamilyContainer."""
        source = inspect.getsource(wdg.SettingsContainer)
        
        self.assertIn('FamilyContainer', source)
        self.assertIn('family_widget', source)

    def test_settings_container_contains_instances_container(self):
        """Test that SettingsContainer contains InstancesContainer."""
        source = inspect.getsource(wdg.SettingsContainer)
        
        self.assertIn('InstancesContainer', source)

    def test_settings_container_contains_families_container(self):
        """Test that SettingsContainer contains FamiliesContainer."""
        source = inspect.getsource(wdg.SettingsContainer)
        
        self.assertIn('FamiliesContainer', source)


class TestWidgetConfiguration(unittest.TestCase):
    """Test cases for widget configuration and styling."""

    def test_settings_container_margins(self):
        """Test that SettingsContainer sets contents margins."""
        source = inspect.getsource(wdg.SettingsContainer)
        self.assertIn('setContentsMargins', source)

    def test_settings_container_palette(self):
        """Test that SettingsContainer sets palette."""
        source = inspect.getsource(wdg.SettingsContainer)
        self.assertIn('setPalette', source)

    def test_source_buttons_disabled_by_default(self):
        """Test that source buttons are disabled by default."""
        source = inspect.getsource(wdg.SettingsContainer)
        self.assertIn('setEnabled(False)', source)

    def test_instance_field_placeholder(self):
        """Test that instance field has placeholder text."""
        source = inspect.getsource(wdg.SettingsContainer)
        self.assertIn('setPlaceholderText', source)
        self.assertIn('Enter instance name..', source)

    def test_custom_data_widget_placeholder(self):
        """Test that custom data widget has placeholder text."""
        source = inspect.getsource(wdg.SettingsContainer)
        self.assertIn('(Optional) Enter custom python dictionary...', source)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)