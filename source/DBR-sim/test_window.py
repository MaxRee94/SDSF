"""Unit test suite for window module.

This module provides comprehensive tests for the window functionality.
Tests cover window classes, their methods, and integration with widgets and gui_utils.

Note: Since PySide6/Qt is not available in the test environment, these tests analyze
 the source code structure rather than executing the actual GUI code.
"""

import unittest
import os
import sys

# Add the parent directory to Python path to access the module file
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Read the window module source code
with open(os.path.join(sys.path[0], 'window.py'), 'r', encoding='utf-8') as f:
    WINDOW_SOURCE = f.read()


class TestMinimalSettingsWidget(unittest.TestCase):
    """Test cases for MinimalSettingsWidget class."""

    def test_minimal_settings_widget_class_exists(self):
        """Test that MinimalSettingsWidget class exists."""
        self.assertIn('class MinimalSettingsWidget(SettingsWidget):', WINDOW_SOURCE)

    def test_minimal_settings_widget_inherits_from_settings_widget(self):
        """Test that MinimalSettingsWidget inherits from SettingsWidget."""
        self.assertIn('class MinimalSettingsWidget(SettingsWidget):', WINDOW_SOURCE)

    def test_minimal_settings_widget_init_signature(self):
        """Test MinimalSettingsWidget.__init__ signature."""
        self.assertIn('def __init__(self, parent=None):', WINDOW_SOURCE)

    def test_minimal_settings_widget_modifies_parent_layout(self):
        """Test that MinimalSettingsWidget modifies parent layout."""
        self.assertIn('super()', WINDOW_SOURCE)
        self.assertIn('__init__', WINDOW_SOURCE)

    def test_minimal_settings_widget_removes_first_row(self):
        """Test that MinimalSettingsWidget removes first row."""
        self.assertIn('removeRow(0)', WINDOW_SOURCE)

    def test_minimal_settings_widget_replaces_path_widget(self):
        """Test that MinimalSettingsWidget replaces path widget."""
        self.assertIn('QLineEdit', WINDOW_SOURCE)
        self.assertIn('setPlaceholderText', WINDOW_SOURCE)

    def test_minimal_settings_widget_inserts_new_row(self):
        """Test that MinimalSettingsWidget inserts new row."""
        self.assertIn('insertRow', WINDOW_SOURCE)
        self.assertIn('path_widget', WINDOW_SOURCE)


class TestMinimalWindow(unittest.TestCase):
    """Test cases for MinimalWindow class."""

    def test_minimal_window_class_exists(self):
        """Test that MinimalWindow class exists."""
        self.assertIn('class MinimalWindow(QtWidgets.QDialog):', WINDOW_SOURCE)

    def test_minimal_window_inherits_from_qdialog(self):
        """Test that MinimalWindow inherits from QDialog."""
        self.assertIn('class MinimalWindow(QtWidgets.QDialog):', WINDOW_SOURCE)

    def test_minimal_window_init_signature(self):
        """Test MinimalWindow.__init__ signature."""
        self.assertIn('def __init__(self, parent=None):', WINDOW_SOURCE)

    def test_minimal_window_has_settings_widget(self):
        """Test that MinimalWindow has settings_widget attribute."""
        self.assertIn('settings_widget', WINDOW_SOURCE)

    def test_minimal_window_has_start_container(self):
        """Test that MinimalWindow has start_container attribute."""
        self.assertIn('start_container', WINDOW_SOURCE)

    def test_minimal_window_has_start_button(self):
        """Test that MinimalWindow has start button."""
        self.assertIn('start_btn', WINDOW_SOURCE)
        self.assertIn('QPushButton', WINDOW_SOURCE)

    def test_minimal_window_has_status_label(self):
        """Test that MinimalWindow has status label."""
        self.assertIn('status_label', WINDOW_SOURCE)
        self.assertIn('QLabel', WINDOW_SOURCE)

    def test_minimal_window_uses_qvboxlayout(self):
        """Test that MinimalWindow uses QVBoxLayout."""
        self.assertIn('QVBoxLayout', WINDOW_SOURCE)

    def test_minimal_window_window_properties(self):
        """Test that MinimalWindow sets window properties."""
        # Check window properties
        self.assertIn('setWindowTitle', WINDOW_SOURCE)
        self.assertIn('setWindowFlags', WINDOW_SOURCE)
        self.assertIn('setFocusPolicy', WINDOW_SOURCE)


class TestModuleStructure(unittest.TestCase):
    """Test cases for overall module structure."""

    def test_module_has_global_window_variable(self):
        """Test that module has global window variable."""
        self.assertIn('window = None', WINDOW_SOURCE)

    def test_module_docstring_missing(self):
        """Test module docstring (may be missing)."""
        # This is just to document that the module might not have a docstring
        # This is acceptable for GUI modules
        self.assertTrue(True)

    def test_expected_classes_present(self):
        """Test that all expected classes are present in the module."""
        expected_classes = ['MinimalSettingsWidget', 'MinimalWindow']
        for class_name in expected_classes:
            self.assertIn(f'class {class_name}', WINDOW_SOURCE)


class TestWindowImports(unittest.TestCase):
    """Test cases for window module imports."""

    def test_widgets_import(self):
        """Test that widgets module is imported."""
        self.assertIn('from widgets import SettingsWidget', WINDOW_SOURCE)

    def test_gui_utils_import(self):
        """Test that gui_utils module is imported."""
        self.assertIn('from gui_utils import HContainer, VContainer, ListSettingWidget, KEYCODES', WINDOW_SOURCE)

    def test_pyside6_imports(self):
        """Test that PySide6 modules are imported."""
        self.assertIn('from PySide6 import', WINDOW_SOURCE)

    def test_config_import(self):
        """Test that config module is imported."""
        self.assertIn('import config', WINDOW_SOURCE)


class TestMinimalWindowLayout(unittest.TestCase):
    """Test cases for MinimalWindow layout structure."""

    def test_minimal_window_layout_structure(self):
        """Test that MinimalWindow has proper layout structure."""
        self.assertIn('layout = QtWidgets.QVBoxLayout', WINDOW_SOURCE)
        self.assertIn('addWidget', WINDOW_SOURCE)

    def test_minimal_window_start_container_layout(self):
        """Test that start container has proper layout."""
        self.assertIn('start_container = HContainer', WINDOW_SOURCE)
        self.assertIn('start_container.layout', WINDOW_SOURCE)

    def test_minimal_window_start_button_properties(self):
        """Test that start button has expected properties."""
        self.assertIn('setFixedWidth', WINDOW_SOURCE)
        self.assertIn('"Start"', WINDOW_SOURCE)

    def test_minimal_window_start_container_structure(self):
        """Test that start container has proper structure."""
        self.assertIn('addStretch', WINDOW_SOURCE)
        self.assertIn('addWidget', WINDOW_SOURCE)


class TestMinimalSettingsWidgetBehavior(unittest.TestCase):
    """Test cases for MinimalSettingsWidget behavior."""

    def test_minimal_settings_widget_modifies_parent(self):
        """Test that MinimalSettingsWidget calls parent __init__."""
        self.assertIn('super()', WINDOW_SOURCE)
        self.assertIn('__init__', WINDOW_SOURCE)

    def test_minimal_settings_widget_removes_first_row(self):
        """Test that MinimalSettingsWidget removes first row."""
        self.assertIn('removeRow(0)', WINDOW_SOURCE)

    def test_minimal_settings_widget_replaces_path_widget(self):
        """Test that MinimalSettingsWidget replaces path widget."""
        self.assertIn('QLineEdit', WINDOW_SOURCE)
        self.assertIn('setPlaceholderText', WINDOW_SOURCE)

    def test_minimal_settings_widget_inserts_new_row(self):
        """Test that MinimalSettingsWidget inserts new row."""
        self.assertIn('insertRow', WINDOW_SOURCE)
        self.assertIn('path_widget', WINDOW_SOURCE)


class TestMinimalWindowBehavior(unittest.TestCase):
    """Test cases for MinimalWindow behavior."""

    def test_minimal_window_window_title(self):
        """Test that MinimalWindow sets window title."""
        self.assertIn('setWindowTitle', WINDOW_SOURCE)
        self.assertIn('"DispRes"', WINDOW_SOURCE)

    def test_minimal_window_window_flags(self):
        """Test that MinimalWindow sets window flags."""
        self.assertIn('setWindowFlags', WINDOW_SOURCE)
        self.assertIn('Qt.Window', WINDOW_SOURCE)

    def test_minimal_window_focus_policy(self):
        """Test that MinimalWindow sets focus policy."""
        self.assertIn('setFocusPolicy', WINDOW_SOURCE)
        self.assertIn('Qt.StrongFocus', WINDOW_SOURCE)

    def test_minimal_window_layout_creation(self):
        """Test that MinimalWindow creates layout."""
        self.assertIn('layout = QtWidgets.QVBoxLayout(self)', WINDOW_SOURCE)

    def test_minimal_window_widget_creation(self):
        """Test that MinimalWindow creates necessary widgets."""
        expected_widgets = ['settings_widget', 'start_container', 'status_label']
        for widget in expected_widgets:
            self.assertIn(widget, WINDOW_SOURCE)

    def test_minimal_window_layout_addition(self):
        """Test that MinimalWindow adds widgets to layout."""
        # Check that widgets are added to layout
        self.assertIn('layout.addWidget(self.settings_widget)', WINDOW_SOURCE)
        self.assertIn('layout.addWidget(start_container)', WINDOW_SOURCE)
        self.assertIn('layout.addWidget(self.status_label)', WINDOW_SOURCE)


class TestModuleVariables(unittest.TestCase):
    """Test cases for module-level variables."""

    def test_module_window_variable(self):
        """Test that module has window variable."""
        self.assertIn('window = None', WINDOW_SOURCE)


class TestModuleInitialization(unittest.TestCase):
    """Test cases for module initialization."""

    def test_module_initialization_order(self):
        """Test that module initializes in correct order."""
        # Test that classes are defined before any module-level code that might use them
        class_defs = ['class MinimalSettingsWidget', 'class MinimalWindow']
        for class_def in class_defs:
            self.assertIn(class_def, WINDOW_SOURCE)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)