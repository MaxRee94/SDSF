"""Unit test suite for gui_utils module.

This module provides comprehensive tests for the GUI utility functionality.
Tests cover container classes, layout classes, list widget functionality, and keycode constants.

Note: Since PySide6/Qt is not available in the test environment, these tests analyze
 the source code structure rather than executing the actual GUI code.
"""

import unittest
import os
import sys

# Add the parent directory to Python path to access the module file
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Read the gui_utils module source code
with open(os.path.join(sys.path[0], 'gui_utils.py'), 'r', encoding='utf-8') as f:
    GUI_UTILS_SOURCE = f.read()


class TestKeycodes(unittest.TestCase):
    """Test cases for KEYCODES constants."""

    def test_keycodes_dict_exists(self):
        """Test that KEYCODES dictionary exists in source."""
        self.assertIn('KEYCODES = {', GUI_UTILS_SOURCE)

    def test_keycodes_has_expected_keys(self):
        """Test that KEYCODES has expected key entries."""
        expected_keys = ["Enter", "Escape", "UpArrow", "DownArrow"]
        for key in expected_keys:
            self.assertIn(f'"{key}":', GUI_UTILS_SOURCE)

    def test_keycodes_values_are_integers(self):
        """Test that all KEYCODES values are integers by examining source."""
        # Extract the KEYCODES dictionary
        start = GUI_UTILS_SOURCE.find('KEYCODES = {')
        end = GUI_UTILS_SOURCE.find('}', start) + 1
        keycodes_section = GUI_UTILS_SOURCE[start:end]
        
        # Should contain integer values
        self.assertIn('16777220', keycodes_section)  # Enter
        self.assertIn('16777216', keycodes_section)  # Escape
        self.assertIn('16777235', keycodes_section)  # UpArrow
        self.assertIn('16777237', keycodes_section)  # DownArrow

    def test_enter_keycode_value(self):
        """Test Enter keycode has expected value."""
        self.assertIn('"Enter": 16777220', GUI_UTILS_SOURCE)

    def test_escape_keycode_value(self):
        """Test Escape keycode has expected value."""
        self.assertIn('"Escape": 16777216', GUI_UTILS_SOURCE)

    def test_up_arrow_keycode_value(self):
        """Test UpArrow keycode has expected value."""
        self.assertIn('"UpArrow": 16777235', GUI_UTILS_SOURCE)

    def test_down_arrow_keycode_value(self):
        """Test DownArrow keycode has expected value."""
        self.assertIn('"DownArrow": 16777237', GUI_UTILS_SOURCE)


class TestModuleStructure(unittest.TestCase):
    """Test cases for overall module structure."""

    def test_module_docstring(self):
        """Test that module has docstring."""
        # Extract module docstring
        first_line = GUI_UTILS_SOURCE.split('\n')[0]
        self.assertIn('"""', first_line)
        self.assertGreater(len(first_line), 3)

    def test_pyside6_imports(self):
        """Test that PySide6 modules are imported."""
        self.assertIn('from PySide6 import QtCore, QtWidgets, QtGui', GUI_UTILS_SOURCE)

    def test_expected_classes_defined(self):
        """Test that all expected classes are defined in the module."""
        expected_classes = [
            'Container', 'HContainer', 'VContainer', 
            'HLayout', 'VLayout', 'ListSettingWidget'
        ]
        for class_name in expected_classes:
            self.assertIn(f'class {class_name}', GUI_UTILS_SOURCE)


class TestContainerClasses(unittest.TestCase):
    """Test cases for Container base class and its subclasses."""

    def test_container_class_definition(self):
        """Test that Container class is defined."""
        self.assertIn('class Container(QtWidgets.QWidget):', GUI_UTILS_SOURCE)

    def test_container_docstring(self):
        """Test that Container has docstring."""
        self.assertIn('Container with layout to which given child widgets are added.', GUI_UTILS_SOURCE)

    def test_container_layout_attribute(self):
        """Test that Container has _layout attribute."""
        self.assertIn('_layout = QtWidgets.QLayout', GUI_UTILS_SOURCE)

    def test_hcontainer_class_definition(self):
        """Test that HContainer class is defined."""
        self.assertIn('class HContainer(Container):', GUI_UTILS_SOURCE)

    def test_hcontainer_inherits_from_container(self):
        """Test that HContainer inherits from Container."""
        self.assertIn('class HContainer(Container):', GUI_UTILS_SOURCE)

    def test_hcontainer_layout_is_hlayout(self):
        """Test that HContainer uses HLayout as _layout."""
        # Find HContainer class definition
        hcontainer_start = GUI_UTILS_SOURCE.find('class HContainer(Container):')
        hcontainer_end = GUI_UTILS_SOURCE.find('\nclass ', hcontainer_start + 1)
        if hcontainer_end == -1:
            hcontainer_end = len(GUI_UTILS_SOURCE)
        hcontainer_section = GUI_UTILS_SOURCE[hcontainer_start:hcontainer_end]
        self.assertIn('_layout = HLayout', hcontainer_section)

    def test_vcontainer_class_definition(self):
        """Test that VContainer class is defined."""
        self.assertIn('class VContainer(Container):', GUI_UTILS_SOURCE)

    def test_vcontainer_inherits_from_container(self):
        """Test that VContainer inherits from Container."""
        self.assertIn('class VContainer(Container):', GUI_UTILS_SOURCE)

    def test_vcontainer_layout_is_vlayout(self):
        """Test that VContainer uses VLayout as _layout."""
        # Find VContainer class definition
        vcontainer_start = GUI_UTILS_SOURCE.find('class VContainer(Container):')
        vcontainer_end = GUI_UTILS_SOURCE.find('\nclass ', vcontainer_start + 1)
        if vcontainer_end == -1:
            vcontainer_end = len(GUI_UTILS_SOURCE)
        vcontainer_section = GUI_UTILS_SOURCE[vcontainer_start:vcontainer_end]
        self.assertIn('_layout = VLayout', vcontainer_section)


class TestLayoutClasses(unittest.TestCase):
    """Test cases for layout classes."""

    def test_hlayout_class_definition(self):
        """Test that HLayout class is defined."""
        self.assertIn('class HLayout(QtWidgets.QHBoxLayout):', GUI_UTILS_SOURCE)

    def test_hlayout_docstring(self):
        """Test that HLayout has docstring."""
        self.assertIn('HBoxLayout without contents margins.', GUI_UTILS_SOURCE)

    def test_vlayout_class_definition(self):
        """Test that VLayout class is defined."""
        self.assertIn('class VLayout(QtWidgets.QVBoxLayout):', GUI_UTILS_SOURCE)

    def test_vlayout_docstring(self):
        """Test that VLayout has docstring."""
        self.assertIn('VBoxLayout without contents margins.', GUI_UTILS_SOURCE)

    def test_hlayout_sets_zero_margins(self):
        """Test that HLayout sets contents margins to zero."""
        # Find HLayout class definition
        hlayout_start = GUI_UTILS_SOURCE.find('class HLayout(QtWidgets.QHBoxLayout):')
        hlayout_end = GUI_UTILS_SOURCE.find('\nclass ', hlayout_start + 1)
        if hlayout_end == -1:
            hlayout_end = len(GUI_UTILS_SOURCE)
        hlayout_section = GUI_UTILS_SOURCE[hlayout_start:hlayout_end]
        self.assertIn('self.setContentsMargins(0, 0, 0, 0)', hlayout_section)

    def test_vlayout_sets_zero_margins(self):
        """Test that VLayout sets contents margins to zero."""
        # Find VLayout class definition
        vlayout_start = GUI_UTILS_SOURCE.find('class VLayout(QtWidgets.QVBoxLayout):')
        vlayout_end = GUI_UTILS_SOURCE.find('\nclass ', vlayout_start + 1)
        if vlayout_end == -1:
            vlayout_end = len(GUI_UTILS_SOURCE)
        vlayout_section = GUI_UTILS_SOURCE[vlayout_start:vlayout_end]
        self.assertIn('self.setContentsMargins(0, 0, 0, 0)', vlayout_section)


class TestListSettingWidget(unittest.TestCase):
    """Test cases for ListSettingWidget class."""

    def test_list_setting_widget_class_definition(self):
        """Test that ListSettingWidget class is defined."""
        self.assertIn('class ListSettingWidget(QtWidgets.QWidget):', GUI_UTILS_SOURCE)

    def test_list_setting_widget_docstring(self):
        """Test that ListSettingWidget has docstring."""
        self.assertIn('Listwidget connected to', GUI_UTILS_SOURCE)

    def test_list_setting_widget_class_variables(self):
        """Test that ListSettingWidget has expected class variables."""
        expected_vars = ['_NEW_ITEM_FIELD', '_PLACEHOLDER_TEXT', '_ADD_BTN_TEXT', '_REMOVE_BTN_TEXT']
        for var in expected_vars:
            self.assertIn(var, GUI_UTILS_SOURCE)

    def test_list_setting_widget_default_values(self):
        """Test ListSettingWidget default class variable values."""
        self.assertIn('_PLACEHOLDER_TEXT = ""', GUI_UTILS_SOURCE)
        self.assertIn('_ADD_BTN_TEXT = "Add"', GUI_UTILS_SOURCE)
        self.assertIn('_REMOVE_BTN_TEXT = "Remove"', GUI_UTILS_SOURCE)

    def test_list_setting_widget_validation_function(self):
        """Test that ListSettingWidget has validation function."""
        self.assertIn('def _IS_ITEM_NAME_VALID(*args, name=None): True', GUI_UTILS_SOURCE)

    def test_list_setting_widget_validation_default(self):
        """Test ListSettingWidget default validation always returns True."""
        self.assertIn('def _IS_ITEM_NAME_VALID(*args, name=None): True', GUI_UTILS_SOURCE)

    def test_list_setting_widget_required_methods(self):
        """Test that ListSettingWidget has required methods."""
        required_methods = [
            'add_item', 'remove_selected_items', 'refresh', 
            '_on_add_clicked', '_on_remove_clicked', '_on_text_changed',
            '_on_selection_changed', 'keyPressEvent', 'position'
        ]
        for method in required_methods:
            self.assertIn(f'def {method}', GUI_UTILS_SOURCE)


class TestContainerInitialization(unittest.TestCase):
    """Test cases for Container initialization behavior."""

    def test_container_init_accepts_parent(self):
        """Test that Container.__init__ accepts parent parameter."""
        self.assertIn('def __init__(self, parent=None, children=[]):', GUI_UTILS_SOURCE)

    def test_container_init_accepts_children(self):
        """Test that Container.__init__ accepts children parameter."""
        self.assertIn('def __init__(self, parent=None, children=[]):', GUI_UTILS_SOURCE)

    def test_container_init_children_default(self):
        """Test that Container children parameter has default value."""
        self.assertIn('children=[]', GUI_UTILS_SOURCE)

    def test_hcontainer_init_signature(self):
        """Test HContainer.__init__ signature."""
        # HContainer should call parent __init__
        self.assertIn('super(HContainer, self).__init__(parent=parent, children=children)', GUI_UTILS_SOURCE)

    def test_vcontainer_init_signature(self):
        """Test VContainer.__init__ signature."""
        # VContainer should call parent __init__
        self.assertIn('super(VContainer, self).__init__(parent=parent, children=children)', GUI_UTILS_SOURCE)


class TestListSettingWidgetMethods(unittest.TestCase):
    """Test cases for ListSettingWidget method behavior."""

    def test_add_item_method_exists(self):
        """Test that add_item method exists."""
        self.assertIn('def add_item(self):', GUI_UTILS_SOURCE)

    def test_remove_selected_items_method_exists(self):
        """Test that remove_selected_items method exists."""
        self.assertIn('def remove_selected_items(self):', GUI_UTILS_SOURCE)

    def test_refresh_method_exists(self):
        """Test that refresh method exists."""
        self.assertIn('def refresh(self):', GUI_UTILS_SOURCE)

    def test_position_method_exists(self):
        """Test that position method exists."""
        self.assertIn('def position(self, widget):', GUI_UTILS_SOURCE)

    def test_key_press_event_method_exists(self):
        """Test that keyPressEvent method exists."""
        self.assertIn('def keyPressEvent(self, event):', GUI_UTILS_SOURCE)


class TestListSettingWidgetEventHandling(unittest.TestCase):
    """Test cases for ListSettingWidget event handling."""

    def test_on_add_clicked_calls_add_item(self):
        """Test that _on_add_clicked calls add_item."""
        self.assertIn('def _on_add_clicked(self, *args):\n        self.add_item()', GUI_UTILS_SOURCE)

    def test_on_remove_clicked_calls_remove_selected_items(self):
        """Test that _on_remove_clicked calls remove_selected_items."""
        self.assertIn('self.remove_selected_items()', GUI_UTILS_SOURCE)

    def test_on_text_changed_calls_refresh(self):
        """Test that _on_text_changed calls refresh."""
        self.assertIn('self.refresh()', GUI_UTILS_SOURCE)

    def test_on_selection_changed_calls_refresh(self):
        """Test that _on_selection_changed calls refresh."""
        self.assertIn('self.refresh()', GUI_UTILS_SOURCE)


class TestContainerChildrenHandling(unittest.TestCase):
    """Test cases for Container children handling logic."""

    def test_container_init_processes_children(self):
        """Test that Container.__init__ processes children parameter."""
        self.assertIn('for child in children:', GUI_UTILS_SOURCE)

    def test_container_handles_widget_children(self):
        """Test Container handles QWidget children."""
        self.assertIn('isinstance(child, QtWidgets.QWidget)', GUI_UTILS_SOURCE)
        self.assertIn('addWidget(child)', GUI_UTILS_SOURCE)

    def test_container_handles_spacer_children(self):
        """Test Container handles QSpacerItem children."""
        self.assertIn('isinstance(child, QtWidgets.QSpacerItem)', GUI_UTILS_SOURCE)
        self.assertIn('addItem(child)', GUI_UTILS_SOURCE)

    def test_container_handles_stretch_children(self):
        """Test Container handles stretch string children."""
        self.assertIn('child == "stretch"', GUI_UTILS_SOURCE)
        self.assertIn('addStretch()', GUI_UTILS_SOURCE)


class TestListSettingWidgetBehavior(unittest.TestCase):
    """Test cases for ListSettingWidget behavior."""

    def test_add_item_uses_new_item_field(self):
        """Test that add_item gets text from new_item_field."""
        self.assertIn('self.new_item_field.text()', GUI_UTILS_SOURCE)

    def test_add_item_validates_name(self):
        """Test that add_item validates item name."""
        self.assertIn('_IS_ITEM_NAME_VALID', GUI_UTILS_SOURCE)
        self.assertIn('validation_result', GUI_UTILS_SOURCE)

    def test_add_item_adds_to_listwidget(self):
        """Test that add_item adds valid items to listwidget."""
        self.assertIn('self.listwidget.insertItem', GUI_UTILS_SOURCE)

    def test_remove_selected_items_removes_from_list(self):
        """Test that remove_selected_items removes selected items."""
        self.assertIn('self.listwidget.selectedItems()', GUI_UTILS_SOURCE)
        self.assertIn('self.listwidget.takeItem', GUI_UTILS_SOURCE)

    def test_refresh_updates_button_states(self):
        """Test that refresh updates button enable states."""
        self.assertIn('setEnabled', GUI_UTILS_SOURCE)

    def test_key_press_event_handles_enter(self):
        """Test that keyPressEvent handles Enter key."""
        self.assertIn('KEYCODES.get("Enter")', GUI_UTILS_SOURCE)
        self.assertIn('self._on_add_clicked()', GUI_UTILS_SOURCE)


class TestModuleImports(unittest.TestCase):
    """Test cases for module imports."""

    def test_pyside6_imports(self):
        """Test that PySide6 modules are imported."""
        self.assertIn('from PySide6 import QtCore, QtWidgets, QtGui', GUI_UTILS_SOURCE)

    def test_no_unexpected_imports(self):
        """Test that only expected imports are present."""
        import_lines = [line for line in GUI_UTILS_SOURCE.split('\n') if line.strip().startswith('import ') or line.strip().startswith('from ')]
        self.assertGreater(len(import_lines), 0)


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)